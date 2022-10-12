import sys
import time
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import progressbar as progbar

from sklearn.linear_model import TheilSenRegressor
from ..pafstats import parse_paf
from ..config import Config
from ..argparse import ArgParser
from ..index import str_to_coord
from ..fast5 import Fast5Reader

from .io import BAM, Guppy

import _uncalled

from ..signal_processor import SignalProcessor

from . import Bcaln, Tracks

import multiprocessing as mp

METHODS = {
    "guided" : "BandedDTW", 
    "static" : "StaticBDTW",
    "global" : "GlobalDTW", 
}


def dtw(conf):
    if conf.tracks.io.guppy_in is None:
        run_dtw(conf)
        return

    guppy_in = Guppy(conf, "r")

    #for config in guppy_in.iter_batches():
    #    print(config.tracks.io.bam_out)

    pbar = progbar.ProgressBar(len(guppy_in.fast5_paths), widgets=[
        progbar.Timer(),
        progbar.Bar(),
        progbar.ETA(),
    ])

    with mp.Pool(processes=conf.threads) as pool:
        sys.stderr.write(f"Using {conf.threads} processes\n")
        pbar.start()
        for i,_ in enumerate(pool.imap_unordered(run_dtw, guppy_in.iter_batches(), chunksize=1)):
            pbar.update(i)

    tracks.close()
    pbar.finish()

def run_dtw(conf):

    """Perform DTW alignment guided by basecalled alignments"""
    conf.fast5_reader.load_bc = True
    conf.export_static()

    tracks = Tracks(conf=conf)

    #tracks.model = tracks.model.get_normalized(*tracks.model.norm_mad_params(tracks.model.means))

    clip_coords = tracks.coords

    fast5s = Fast5Reader(conf=conf)

    read_filter = fast5s.get_read_filter()

    if not isinstance(tracks.input, BAM):
        raise ValueError("Must specify BAM input")

    sigproc = SignalProcessor(tracks[tracks.output_track].model, conf)

    n_reads = 0

    for read in fast5s:
        aligned = False
        for aln in tracks.input.get_alns(read.id):
            #sys.stderr.write(f"{read.id}\n")

            dtw = GuidedDTW(tracks, sigproc, read, aln, conf)

            if dtw.df is None:
                sys.stderr.write(f"Warning: {read.id} failed\n")
                continue

            if conf.bc_cmp:
                tracks.calc_compare("bcaln", True, True)

            tracks.write_alignment()

            aligned = True

        if aligned:
            #pbar.update(n_reads)
            n_reads += 1

    tracks.close()

    #pbar.finish()

class GuidedDTW:

    #TODO do more in constructor using prms, not in main
    #def __init__(self, track, read, paf, conf=None, **kwargs):
    def __init__(self, tracks, sigproc, read, aln, conf=None, **kwargs):
        self.conf = read.conf if conf is None else conf
        self.prms = self.conf.dtw

        bcaln = Bcaln(conf, tracks.index, read, aln, tracks.coords)
        if bcaln.empty:
            self.df = None
            return
        
        signal = sigproc.process(read, False)

        #if self.conf.normalizer.mode == "ref_mom":
        #    signal.normalize(

        aln_id, self.coords = tracks.init_alignment(read.id, read.filename, bcaln.coords, {"bcaln" : bcaln.df}, read=signal)

        self.index = tracks.index
        self.model = sigproc.model

        self.bcaln = bcaln.df#[bcaln.df["indel"] == 0]

        self.ref_gaps = list(sorted(bcaln.ref_gaps))

        kmer_blocks = list()
        gap_lens = list()
        ref_shift = 0

        mref_min = self.coords.mrefs.min()#-self.index.trim[0]
        mref_max = self.coords.mrefs.max()+1#+self.index.trim[1]

        self.block_coords = list()
        block_st = mref_min
        shift = block_st#+5
        for gap_st, gap_en in self.ref_gaps:
            self.block_coords.append([block_st, gap_st, shift])
            block_st = gap_en
            shift += gap_en-gap_st
        self.block_coords.append([block_st,mref_max,shift])

        #kmer_blocks = [kmers.loc[st:en] for st,en,_ in self.block_coords]
        new_kmers = self.index.get_kmers([(s,e) for s,e,_ in self.block_coords], conf.is_rna)

        self.block_coords[0][0] += self.index.trim[1]
        self.block_coords[-1][1] -= self.index.trim[0]
        mrefs = np.concatenate([np.arange(s,e) for s,e,_ in self.block_coords])

        #self.ref_kmers = pd.concat(kmer_blocks)
        self.ref_kmers = pd.Series(new_kmers, index=mrefs)

        ref_means = self.model[self.ref_kmers]

        self.method = self.prms.band_mode
        if not self.method in METHODS:
            opts = "\", \"".join(METHODS.keys())
            raise ValueError(f"Error: unrecongized DTW method \"{method}. Must be one of \"{opts}\"")

        method = METHODS[self.method]
        self.dtw_fn = getattr(_uncalled, f"{method}K{self.model.K}", None)

        if self.dtw_fn is None:
            raise ValueError(f"Invalid DTW k-mer length {self.model.K}")

        self.samp_min = self.bcaln["start"].min()
        self.samp_max = self.bcaln["start"].max()
        self.evt_start, self.evt_end = signal.event_bounds(self.samp_min, self.samp_max)

        if self.conf.normalizer.mode == "ref_mom":
            
            if self.conf.normalizer.median:
                med = np.median(ref_means)
                #tgt = (med, np.median(np.abs(med-ref_means)))
                tgt = (med, np.median(np.abs(med-ref_means)))
            else:
                tgt = (ref_means.mean(), ref_means.std())


            signal.normalize_mom(ref_means.mean(), ref_means.std())#, self.evt_start, self.evt_end)
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

        df = self._calc_dtw(signal)

        for i in range(self.prms.iterations-1):
            reg = self.renormalize(signal, df)
            signal.normalize(reg.coef_, reg.intercept_)
            df = self._calc_dtw(signal)


        df["kmer"] = self.ref_kmers.loc[df["mref"]].to_numpy()
        self.df = df.set_index("mref")

        tracks.write_dtw_events(self.df, aln_id=aln_id, read=signal)

        if self.bands is not None:
            tracks.add_layers("band", self.bands, aln_id=aln_id)

        self.empty = False

    def renormalize(self, signal, aln):
        kmers = self.ref_kmers[aln["mref"]]
        model_current = self.model[kmers]
        reg = TheilSenRegressor(random_state=0)
        return reg.fit(aln[["current"]], model_current)

    #TODO refactor inner loop to call function per-block
    def _calc_dtw(self, signal):
        self.mats = list()

        path_qrys = list()
        path_refs = list()

        #mref_st = self.ref_kmers.index[0]
        mref_st = self.bcaln.index[self.bcaln.index.searchsorted(self.ref_kmers.index[0])]
        mref_en = self.bcaln.index[self.bcaln.index.searchsorted(self.ref_kmers.index[-1])]
        #mref_en = self.bcaln.index[-1]+1

        #kmers = self.ref_kmers.loc[mref_st:mref_en]
        #kmers = self.ref_kmers

        samp_st = self.bcaln.loc[mref_st]["start"]
        samp_en = self.bcaln.loc[mref_en]["start"] + self.bcaln.loc[mref_en]["length"]
        #samp_st = self.bcaln.iloc[0]["start"]
        #samp_en = self.bcaln.iloc[-1]["start"] + self.bcaln.iloc[-1]["length"]

        read_block = signal.to_df()[self.evt_start:self.evt_end]#.sample_range(samp_st, samp_en)

        block_signal = read_block['mean'].to_numpy()
        
        args = self._get_dtw_args(read_block, self.ref_kmers)


        dtw = self.dtw_fn(*args)

        
        path = np.flip(dtw.path)
        #TODO shouldn't need to clip, error in bdtw
        evts = read_block.index[np.clip(path['qry'], 0, len(read_block))]
        mrefs = self.ref_kmers.index[path['ref']]

        if self.prms.save_bands and hasattr(dtw, "ll"):
            self.bands = self._ll_to_df(dtw.ll, read_block, mrefs[0], len(self.ref_kmers))
        else:
            self.bands = None

        return pd.DataFrame({'mref': mrefs}, 
                            index = evts,
                            dtype='Int64') \
                 .join(signal.to_df()) \
                 .drop(columns=['mask'], errors='ignore') \
                 .rename(columns={'mean' : 'current', 'stdv' : 'current_stdv'})

    def _get_dtw_args(self, read_block, ref_kmers):
        qry_len = len(read_block)
        ref_len = len(ref_kmers)

        #should maybe move to C++
        if self.method == "guided":
            band_count = qry_len + ref_len

            shift = int(np.round(self.prms.band_shift*self.prms.band_width))

            ar = _uncalled.PyArrayI32
            #mrefs = np.array(ref_kmers.index) #np.array(self.bcaln.index)
            aln = self.bcaln[self.bcaln.index.isin(ref_kmers.index)]
            mrefs = np.array(aln.index)
            for st,en,sh in self.block_coords:
                mrefs[aln.index.isin(pd.RangeIndex(st,en))] -= sh

            bands = _uncalled.get_guided_bands(ar(mrefs), ar(aln["start"]), ar(read_block['start']), band_count, shift)
            
            return (self.prms, _uncalled.PyArrayF32(read_block['mean']), self.model.kmer_array(ref_kmers), self.model.instance, _uncalled.PyArrayCoord(bands))

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
