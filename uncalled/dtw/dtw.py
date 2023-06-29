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
from scipy.stats import linregress
from ..config import Config
from ..argparse import ArgParser
from ..index import str_to_coord
from .. import ExceptionWrapper

import _uncalled
from _uncalled import EventDetector

from ..signal_processor import ProcessedRead

from . import Tracks
from ..aln import AlnDF

from _uncalled import ReadBuffer

import multiprocessing as mp


import memray


#from concurrent.futures import ProcessPoolExecutor as Pool

#from https://stackoverflow.com/questions/6126007/python-getting-a-traceback-from-a-multiprocessing-process


METHODS = {
    "guided" : "BandedDTW", 
    "static" : "StaticBDTW",
    "global" : "GlobalDTW", 
}


def dtw(conf):
    conf.read_index.load_signal = True
    conf.tracks.layers.append("moves")
    #conf.export_static()

    #if len(conf.read_index.read_index) == 0:
    #    raise ValueError("Must specify fast5 index (-x/--fast5-index)")

    if conf.tracks.io.processes == 1:
        dtw_single(conf)
    else:
        dtw_pool(conf)

def dtw_pool(conf):
    #mp.set_start_method("spawn")
    t = time()
    tracks = Tracks(conf=conf)
    assert(tracks.output is not None)
    #tracks.output.set_model(tracks.model)
    i = 0
    _ = tracks.read_index.default_model #load property
    for chunk in dtw_pool_iter(tracks):
        i += len(chunk)
        tracks.output.write_buffer(chunk)

def dtw_pool_iter(tracks):
    def iter_args(): 
        i = 0
        aln_count = 0
        for read_ids, bams in tracks.bam_in.iter_str_chunks():
            reads = tracks.read_index.subset(read_ids)
            yield (tracks.conf, tracks.model, bams, reads, aln_count, tracks.bam_in.header)
            aln_count += len(bams)
    try:
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=tracks.conf.tracks.io.processes, maxtasksperchild=4) as pool:


            i = 0
            if tracks.conf.tracks.io.ordered_out:
                itr = pool.imap
            else:
                itr = pool.imap_unordered

            for out in itr(dtw_worker, iter_args(), chunksize=1):
                i += len(out)
                yield out 

    except Exception as e:
        raise ExceptionWrapper(e).re_raise()

EVDT_PRESETS = {
    450.0 : EventDetector.PRMS_450BPS,
    70.0 : EventDetector.PRMS_70BPS,
}

def init_model(tracks):
    if tracks.model is None:
        raise RuntimeError("Failed to auto-detect pore model. Please specify --pore-model flag (and add --rna if aligning RNA reads)")

    evdt = EVDT_PRESETS.get(tracks.model.bases_per_sec, None)
    tracks.conf.load_group("event_detector", evdt, keep_nondefaults=True)


def dtw_worker(p):
    conf,model,bams,reads,aln_start,header = p

    conf.tracks.io.buffered = True
    #conf.tracks.io.bam_in = None

    header = pysam.AlignmentHeader.from_dict(header)

    tracks = Tracks(model=model, read_index=reads, conf=conf)
    init_model(tracks)

    #tracks.norm_params = pd.read_csv("/scratch1/skovaka/curlcakes/unm/nanopolish/eventalign/npl_unm_cov100_smry.txt", sep="\t").set_index("read_name").sort_index()
    #tracks.norm_params = pd.read_csv("/scratch1/skovaka/curlcakes/m6a/nanopolish/eventalign/npl_m6a_cov100_smry.txt", sep="\t").set_index("read_name").sort_index()

    i = 0
    for bam in bams:
        bam = pysam.AlignedSegment.fromstring(bam, header)
        aln = tracks.bam_in.sam_to_aln(bam)
        aln.instance.id += aln_start
        if aln is not None:
            dtw = GuidedDTW(tracks, aln)
        else:
            sys.stderr.write(f"Warning: {aln.read_id} BAM parse failed\n")

        i += 1

    tracks.close()

    return tracks.output.out_buffer

def dtw_single(conf):
    """Perform DTW alignment guided by basecalled alignments"""

    tracks = Tracks(conf=conf)
    init_model(tracks)
    sys.stderr.write(f"Using model '{tracks.model.name}'\n")

    assert(tracks.output is not None)

    for aln in tracks.bam_in.iter_alns():
        dtw = GuidedDTW(tracks, aln)

    tracks.close()

class GuidedDTW:

    def process(self, read):
        evdt = EventDetector(self.conf.event_detector)
        return ProcessedRead(evdt.process_read(read))

    #def __init__(self, tracks, bam, **kwargs):
    def __init__(self, tracks, aln):
        self.conf = tracks.conf
        self.prms = self.conf.dtw

        self.aln = aln
        self.moves = self.aln.moves #moves.aln
        #print(self.moves.to_pandas())
        if self.moves.empty():
            sys.stderr.write(f"Warning: moves missing for read {aln.read_id}\n")
            return

        read = self.aln.read
        #read = ReadBuffer(read.id,read.channel,read.number,read.start,read.signal)


        signal = self.process(read)

        self.model = tracks.model

        self.seq = self.aln.seq
        #print(len(self.aln.moves))
        #print(len(self.seq))

        self.method = self.prms.band_mode
        if not self.method in METHODS:
            opts = "\", \"".join(METHODS.keys())
            raise ValueError(f"Error: unrecongized DTW method \"{method}. Must be one of \"{opts}\"")

        method = METHODS[self.method]

        self.dtw_fn = None
        if isinstance(self.model.instance, _uncalled.PoreModelU16):
            self.dtw_fn = getattr(_uncalled, f"{method}U16", None)
        elif isinstance(self.model.instance, _uncalled.PoreModelU32):
            self.dtw_fn = getattr(_uncalled, f"{method}U32", None)
        if self.dtw_fn is None:
            raise ValueError(f"Unknown PoreModel type: {model.instance}")

        self.samp_min = self.moves.samples.starts[0]
        self.samp_max = self.moves.samples.ends[len(self.moves)-1]
        self.evt_start, self.evt_end = signal.event_bounds(self.samp_min, self.samp_max)

        if self.prms.iterations == 0:
            tgt = (0, 1)
        elif self.conf.normalizer.mode == "ref_mom":

            ref_means = self.seq.current.to_numpy()  #self.model[self.ref_kmers]
            
            if self.conf.normalizer.median:
                med = np.median(ref_means)
                tgt = (med, np.median(np.abs(med-ref_means)))
            else:
                tgt = (ref_means.mean(), ref_means.std())

        elif self.conf.normalizer.mode == "model_mom":
            if self.conf.normalizer.median:
                med = self.model.means.median()
                tgt = (med, np.median(np.abs(med-self.model.means)))
            else:
                tgt = (self.model.model_mean, self.model.model_stdv)
        else:
            raise ValueError(f"Unknown normalization mode: {self.prms.norm_mode}")

        #scale = 1/self.aln.sam.get_tag("sd")
        #shift = -self.aln.sam.get_tag("sm") * scale
        #signal.normalize(scale, shift)

        #if not aln.read_id in tracks.norm_params.index:
        #    return
        #scale, shift = tracks.norm_params.loc[aln.read_id, ["scale","shift"]]
        #signal.normalize(scale, shift/scale)
        #signal.normalize(scale, shift)

        if self.conf.normalizer.full_read:
            signal.normalize_mom(*tgt)
        else:
            signal.normalize_mom(*tgt, self.evt_start, self.evt_end)

        tracks.set_read(signal)

        if self.prms.iterations > 0:
            success = self._calc_dtw(signal)

            for i in range(self.prms.iterations-1):
                if self.renormalize(signal, self.aln):
                    success = self._calc_dtw(signal)
                else:
                    success = False
                    break

        else:
            starts = self.moves.samples.starts.to_numpy()
            lengths = self.moves.samples.lengths.to_numpy()
            #print(self.aln.read_id)
            dtw = AlnDF(self.seq, starts, lengths)
            dtw.set_signal(signal)
            self.aln.set_dtw(dtw)
            success = True

        
        if success:
            #print("ERE")
            #try:
            #if self.conf.mvcmp and self.aln.mvcmp.empty():
            if self.aln.mvcmp.empty():
                self.aln.calc_mvcmp()

            mask = self.aln.dtw.na_mask

            #mask &= np.array(self.aln.mvcmp.dist) < 5

            if (not self.prms.unmask_splice) and aln.seq.is_spliced():
                mask &= np.array(aln.seq.splice_mask)

            if np.sum(mask) < self.conf.tracks.min_aln_length:
                return

            mask[0] = mask[-1] = True

            self.aln.dtw.mask(mask)

            tracks.write_alignment(self.aln)
            self.empty = False
            return
            #except:
            #    pass
        #else:
        #    print("FAIL")
        sys.stderr.write(f"Failed to write alignment for {read.id}\n")


    def renormalize(self, signal, aln):
        #kmers = self.aln.seq.kmer #self.ref_kmers[aln["mpos"]]
        #model_current = self.model[kmers]
        if aln.mvcmp.empty():
            aln.calc_mvcmp()

        na_mask = np.array(aln.dtw.na_mask)
        dist_mask = np.array(aln.mvcmp.dist) <= self.conf.tracks.min_norm_dist

        #print(np.array(aln.mvcmp.dist))

        mask = na_mask & dist_mask
        if np.sum(mask) < self.conf.tracks.min_aln_length:
            if np.sum(na_mask) < self.conf.tracks.min_aln_length:
                return False
            #sys.stderr.write(f"{aln.read_id}\tnofilter\t{np.sum(mask)}\t{len(mask)}\n")
            mask = na_mask
        
        ref_current = aln.seq.current.to_numpy()[mask] #self.ref_kmers[aln["mpos"]]
        read_current = aln.dtw.current.to_numpy()[mask]
        lr = linregress(read_current, ref_current)
        #reg = TheilSenRegressor(random_state=0)
        #fit = reg.fit(read_current.reshape((-1,1)), ref_current)
        scale, shift = (lr.slope, lr.intercept)
        signal.normalize(scale, shift)
        return True

        #return(fit.coef_[0], fit.intercept_)

    #TODO refactor inner loop to call function per-block
    def _calc_dtw(self, signal):
        read_block = signal.events[self.evt_start:self.evt_end] #.to_df()[self.evt_start:self.evt_end]

        qry_len = self.evt_end - self.evt_start
        ref_len = len(self.seq)

        band_count = qry_len + ref_len
        shift = int(np.round(self.prms.band_shift*self.prms.band_width))
        mv_starts = self.moves.samples.starts.to_numpy()

        bands = _uncalled.get_guided_bands(np.arange(len(self.seq)), mv_starts, read_block['start'], band_count, shift)


        dtw = self.dtw_fn(self.prms, signal, self.evt_start, self.evt_end, self.model.kmer_array(self.seq.kmer.to_numpy()), self.model.instance, _uncalled.PyArrayCoord(bands))

        if np.any(dtw.path["qry"] < self.evt_start) or np.any(dtw.path["ref"] < 0):
            #print("FAAAIIL", dtw.path)
            return False

        dtw.fill_aln(self.aln.instance, self.conf.tracks.count_events)

        if self.conf.tracks.mask_skips is not None:
            self.aln.mask_skips(self.conf.tracks.mask_skips == "keep_best")
        
        #print("HRERE")

        return True


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
