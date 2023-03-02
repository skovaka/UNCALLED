import sys
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import progressbar as progbar
import pysam
from time import time
import os

from sklearn.linear_model import TheilSenRegressor
from ..pafstats import parse_paf
from ..config import Config
from ..argparse import ArgParser
from ..index import str_to_coord

from .io import BAM, Guppy

import _uncalled

from ..signal_processor import ProcessedRead

from . import Bcaln, Tracks

import multiprocessing as mp
#from concurrent.futures import ProcessPoolExecutor as Pool

METHODS = {
    "guided" : "BandedDTW", 
    "static" : "StaticBDTW",
    "global" : "GlobalDTW", 
}


def dtw(conf):
    conf.fast5_reader.load_bc = True
    conf.tracks.load_fast5s = True
    conf.export_static()

    #if len(conf.fast5_reader.fast5_index) == 0:
    #    raise ValueError("Must specify fast5 index (-x/--fast5-index)")

    if conf.tracks.io.processes == 1:
        dtw_single(conf)
    else:
        dtw_pool(conf)

def dtw_pool(conf):
    tracks = Tracks(conf=conf)
    assert(tracks.output is not None)
    #tracks.output.set_model(tracks.model)
    i = 0
    for chunk in dtw_pool_iter(tracks):
        i += len(chunk)
        tracks.output.write_buffer(chunk)
    tracks.close()

def dtw_pool_iter(tracks):
    def iter_args(): 
        i = 0
        for read_ids, bams in tracks.bam_in.iter_str_chunks():
            reads = tracks.read_index.subset(read_ids)
            i += len(bams)
            yield (tracks.conf, bams, reads, tracks.bam_in.header)

    with mp.Pool(processes=tracks.conf.tracks.io.processes) as pool:
        i = 0
        for out in pool.imap(dtw_worker, iter_args(), chunksize=1):
        #for out in pool.imap(dtw_worker, iter_args(), chunksize=1):
            i += len(out)
            yield out 

def dtw_worker(p):
    conf,bams,reads,header = p

    conf.tracks.io.buffered = True
    conf.tracks.io.bam_in = None

    header = pysam.AlignmentHeader.from_dict(header)

    tracks = Tracks(read_index=reads, conf=conf)

    i = 0
    for bam in bams:
        bam = pysam.AlignedSegment.fromstring(bam, header)
        dtw = GuidedDTW(tracks, bam)
        i += 1

    tracks.close()

    return tracks.output.out_buffer


def dtw_single(conf):
    """Perform DTW alignment guided by basecalled alignments"""

    tracks = Tracks(conf=conf)

    assert(tracks.output is not None)

    #for bam in tracks.bam_in.iter_sam():
    #    dtw = GuidedDTW(tracks, bam)

    for sam, aln in tracks.bam_in.iter_moves():
        dtw = GuidedDTW(tracks, aln, sam)

    tracks.close()

class GuidedDTW:

    def process(self, read):
        evdt = _uncalled.EventDetector(self.conf.event_detector)
        return ProcessedRead(evdt.process_read(read))

    #def __init__(self, tracks, bam, **kwargs):
    def __init__(self, tracks, aln, sam):
        self.conf = tracks.conf
        self.prms = self.conf.dtw

        self.aln = aln

        try:
            #read = tracks.read_index[bam.query_name]
            read = tracks.read_index[self.aln.read_id] #TODO don't load twice (reuse from BAM IO)
        except:
            sys.stderr.write(f"Warning: failed to load signal from read {bam.query_name}\n")
            return

        signal = self.process(read)

        self.index = tracks.index
        self.model = tracks.model

        self.moves = self.aln.moves #moves.aln

        #TODO change to tracks.init_alignment(read_id, ref_id, coords)
        #returns python Alignment object (wrapper for AlignmentK*)
        #also store alignment in tracks output buffer

        self.seq = tracks.get_seq(sam.reference_id, self.moves.index)

        self.method = self.prms.band_mode
        if not self.method in METHODS:
            opts = "\", \"".join(METHODS.keys())
            raise ValueError(f"Error: unrecongized DTW method \"{method}. Must be one of \"{opts}\"")

        method = METHODS[self.method]
        self.dtw_fn = getattr(_uncalled, f"{method}K{self.model.K}", None)

        if self.dtw_fn is None:
            raise ValueError(f"Invalid DTW k-mer length {self.model.K}")

        self.samp_min = self.moves.samples.starts[0]
        self.samp_max = self.moves.samples.ends[len(self.moves)-1]
        self.evt_start, self.evt_end = signal.event_bounds(self.samp_min, self.samp_max)

        if self.conf.normalizer.mode == "ref_mom":

            ref_means = self.seq.current.to_numpy()  #self.model[self.ref_kmers]
            
            if self.conf.normalizer.median:
                med = np.median(ref_means)
                #tgt = (med, np.median(np.abs(med-ref_means)))
                tgt = (med, np.median(np.abs(med-ref_means)))
            else:
                tgt = (ref_means.mean(), ref_means.std())

            #signal.normalize_mom(ref_means.mean(), ref_means.std())#, self.evt_start, self.evt_end)
        elif self.conf.normalizer.mode == "model_mom":
            if self.conf.normalizer.median:
                med = self.model.means.median()
                tgt = (med, np.median(np.abs(med-self.model.means)))
            else:
                tgt = (self.model.model_mean, self.model.model_stdv)
        else:
            raise ValueError(f"Unknown normalization mode: {self.prms.norm_mode}")

        if self.conf.normalizer.full_read:
            signal.normalize_mom(*tgt)
        else:
            signal.normalize_mom(*tgt, self.evt_start, self.evt_end)
        
        tracks.set_read(signal)

        if self.prms.iterations > 0:
            self._calc_dtw(signal)

            for i in range(self.prms.iterations-1):
                reg = self.renormalize(signal, df)
                signal.normalize(reg.coef_, reg.intercept_)
                self._calc_dtw(signal)
        else:
            starts = self.moves.samples.starts.to_numpy()
            lengths = self.moves.samples.lengths.to_numpy()
            cur = std = np.array([])
            df = DtwDF(starts, lengths, cur, std)
            df.set_signal(signal)
            df = df.to_df()
            df["mref"] = self.moves.index.expand()

        #if df is None:
        #    self.empty = True
        #    sys.stderr.write(f"Warning: dtw failed for read {read.id}\n")
        #    return
        #tracks.write_dtw_events(self.df, read=signal)#, aln_id=aln_id

        #if self.bands is not None:
        #    tracks.add_layers("band", self.bands)#, aln_id=aln_id)

        #if self.conf.mvcmp:
        #    self.aln.calc_mvcmp()
            #tracks.calc_compare("moves", True, True)

        tracks.write_alignment(self.aln)
        #aln_id, aln_coords = self.init_alignment(aln.read_id, read.filename, moves.coords, {"moves" : moves.df}, bam=bam) #, read=signal

        self.empty = False

    def renormalize(self, signal, aln):
        kmers = self.seq.get_kmer(aln["mref"]) #self.ref_kmers[aln["mref"]]
        model_current = self.model[kmers]
        reg = TheilSenRegressor(random_state=0)
        return reg.fit(aln[["current"]], model_current)

    #TODO refactor inner loop to call function per-block
    def _calc_dtw(self, signal):
        read_block = signal.events[self.evt_start:self.evt_end] #.to_df()[self.evt_start:self.evt_end]

        #print(read_block)
        
        #prms, means, kmers, inst, bands = self._get_dtw_args(read_block)#, self.ref_kmers)

        qry_len = self.evt_end - self.evt_start
        ref_len = len(self.seq)

        band_count = qry_len + ref_len
        shift = int(np.round(self.prms.band_shift*self.prms.band_width))
        mv_starts = self.moves.samples.starts.to_numpy()

        bands = _uncalled.get_guided_bands(np.arange(len(self.seq)), mv_starts, read_block['start'], band_count, shift)

        #return (self.prms, _uncalled.PyArrayF32(read_block['mean']), self.model.kmer_array(self.seq.kmer.to_numpy()), self.model.instance, _uncalled.PyArrayCoord(bands))


        dtw = self.dtw_fn(self.prms, signal, self.evt_start, self.evt_end, self.model.kmer_array(self.seq.kmer.to_numpy()), self.model.instance, _uncalled.PyArrayCoord(bands))

        if np.any(dtw.path["qry"] < 0) or np.any(dtw.path["ref"] < 0):
            return None

        #print(list(dtw.path))

        dtw.fill_aln(self.aln.instance)

        #print(self.aln.dtw.to_pandas().to_string())

    def _get_dtw_args(self, read_block):#, ref_kmers):
        qry_len = len(read_block)
        ref_len = len(self.seq)

        #should maybe move to C++
        if self.method == "guided":
            band_count = qry_len + ref_len

            shift = int(np.round(self.prms.band_shift*self.prms.band_width))
            mv_starts = self.moves.samples.starts.to_numpy()
            bands = _uncalled.get_guided_bands(np.arange(len(self.seq)), mv_starts, read_block['start'], band_count, shift)

            return (self.prms, _uncalled.PyArrayF32(read_block['mean']), self.model.kmer_array(self.seq.kmer.to_numpy()), self.model.instance, _uncalled.PyArrayCoord(bands))

        #elif self.method == "static":
        #    return common

        else:
            return (self.prms, read_block['mean'].to_numpy(), self.model.kmer_array(ref_kmers), self.model.instance, DTW_GLOB)

    #TODO store in ReadAln metadata
    def _ll_to_df(self, ll, read_block, min_mref, ref_len):
        qry_st = np.clip(ll['qry'],                 0, len(read_block)-1)
        qry_en = np.clip(ll['qry']-self.prms.band_width, 0, len(read_block)-1)
        mref_st = min_mref + np.clip(ll['ref'],                 0, ref_len-1)
        mref_en = min_mref + np.clip(ll['ref']+self.prms.band_width, 0, ref_len-1)

        pac_st = self.index.mref_to_pac(mref_en)
        pac_en = self.index.mref_to_pac(mref_st)

        sample_starts = read_block['start'].iloc[qry_en].to_numpy()
        sample_ends = read_block["start"].iloc[qry_st].to_numpy() + read_block["length"].iloc[qry_st].to_numpy()

        grp = pd.DataFrame({
                "pac" : pac_st,
                "pac_end" : pac_en,
                "sample_start" : sample_starts,
                "sample_end" : sample_ends,
              }).groupby("pac")

        df = pd.DataFrame({
            "pac_end" : grp["pac_end"].first(),
            "sample_start" : grp["sample_start"].max(),
            "sample_end" : grp["sample_end"].min(),
        })

        return df

        #return pd.DataFrame({
        #    'samp': band_samps, 
        #    'ref_st': ref_st + band_ref_st, 
        #    'ref_en': ref_st + band_ref_en
        #})

#def collapse_events(dtw, kmer_str=False, start_col="start", length_col="length", mean_col="current", kmer_col="kmer", mask_skips=False):
#
#    dtw["cuml_mean"] = dtw[length_col] * dtw[mean_col]
#
#    grp = dtw.groupby("mref")
#
#    if kmer_col in dtw:
#        if kmer_str:
#            kmers = [nt.kmer_rev(nt.str_to_kmer(k,0)) for k in grp[kmer_col].first()]
#        else:
#            kmers = grp[kmer_col].first()
#    else:
#        kmers = None
#
#    mrefs = grp["mref"].first()
#
#    lengths = grp[length_col].sum()
#
#    dtw = pd.DataFrame({
#        "mref"    : mrefs.astype("int64"),
#        "start"  : grp[start_col].min().astype("uint32"),
#        "length" : lengths.astype("uint32"),
#        "current"   : grp["cuml_mean"].sum() / lengths
#    })
#
#    if kmers is not None:
#        dtw["kmer"] = kmers.astype("uint16")
#
#    if mask_skips:
#        dtw = dtw[~dtw.duplicated("start", False)]
#
#    return dtw.set_index("mref").sort_index()
