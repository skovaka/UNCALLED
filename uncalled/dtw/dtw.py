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
from ..config import Config
from ..argparse import ArgParser
from ..index import str_to_coord

import _uncalled

from ..signal_processor import ProcessedRead

from . import Tracks
from ..aln import AlnDF

import multiprocessing as mp
#from concurrent.futures import ProcessPoolExecutor as Pool

#from https://stackoverflow.com/questions/6126007/python-getting-a-traceback-from-a-multiprocessing-process
import tblib.pickling_support
tblib.pickling_support.install()
class ExceptionWrapper(object):
    def __init__(self, ee):
        self.ee = ee
        __, __, self.tb = sys.exc_info()
    def re_raise(self):
        raise self.ee.with_traceback(self.tb)


METHODS = {
    "guided" : "BandedDTW", 
    "static" : "StaticBDTW",
    "global" : "GlobalDTW", 
}


def dtw(conf):
    conf.tracks.load_fast5s = True
    #conf.export_static()

    #if len(conf.read_index.read_index) == 0:
    #    raise ValueError("Must specify fast5 index (-x/--fast5-index)")

    if conf.tracks.io.processes == 1:
        dtw_single(conf)
    else:
        dtw_pool(conf)

def dtw_pool(conf):
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
        with mp.Pool(processes=tracks.conf.tracks.io.processes) as pool:
            i = 0
            for out in pool.imap(dtw_worker, iter_args(), chunksize=1):
            #for out in pool.imap(dtw_worker, iter_args(), chunksize=1):
                i += len(out)
                yield out 
    except Exception as e:
        raise ExceptionWrapper(e).re_raise()

def dtw_worker(p):
    conf,model,bams,reads,aln_start,header = p

    conf.tracks.io.buffered = True
    #conf.tracks.io.bam_in = None

    header = pysam.AlignmentHeader.from_dict(header)

    tracks = Tracks(model=model, read_index=reads, conf=conf)

    i = 0
    for bam in bams:
        bam = pysam.AlignedSegment.fromstring(bam, header)
        aln = tracks.bam_in.sam_to_aln(bam)
        #sys.stderr.write(f"a {aln.id}\n")
        aln.id += aln_start
        #sys.stderr.write(f"b {aln.id}\n")
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

    assert(tracks.output is not None)

    #for bam in tracks.bam_in.iter_sam():
    #    dtw = GuidedDTW(tracks, bam)

    for aln in tracks.bam_in.iter_alns():
        dtw = GuidedDTW(tracks, aln)

    tracks.close()

class GuidedDTW:

    def process(self, read):
        evdt = _uncalled.EventDetector(self.conf.event_detector)
        return ProcessedRead(evdt.process_read(read))

    #def __init__(self, tracks, bam, **kwargs):
    def __init__(self, tracks, aln):
        self.conf = tracks.conf
        self.prms = self.conf.dtw

        self.aln = aln

        read = self.aln.read

        signal = self.process(read)

        self.model = tracks.model

        self.moves = self.aln.moves #moves.aln
        self.seq = self.aln.seq

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
            success = self._calc_dtw(signal)

            for i in range(self.prms.iterations-1):
                reg = self.renormalize(signal, df)
                signal.normalize(reg.coef_, reg.intercept_)
                success = self._calc_dtw(signal)
        else:
            starts = self.moves.samples.starts.to_numpy()
            lengths = self.moves.samples.lengths.to_numpy()
            cur = std = np.array([])
            dtw = AlnDF(self.seq, starts, lengths)
            dtw.set_signal(signal)
            self.aln.set_dtw(dtw)
            success = True

        if success:
            #try:
            if self.conf.mvcmp:
                self.aln.calc_mvcmp()
            tracks.write_alignment(self.aln)
            self.empty = False
            return
            #except:
            #    pass
        sys.stderr.write(f"Failed to write alignment for {read.id}\n")


    def renormalize(self, signal, aln):
        kmers = self.seq.get_kmer(aln["mpos"]) #self.ref_kmers[aln["mpos"]]
        model_current = self.model[kmers]
        reg = TheilSenRegressor(random_state=0)
        return reg.fit(aln[["current"]], model_current)

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

        if np.any(dtw.path["qry"] < 0) or np.any(dtw.path["ref"] < 0):
            return False

        dtw.fill_aln(self.aln.instance)
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
