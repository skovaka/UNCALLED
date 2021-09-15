import sys
import time
import re
import numpy as np
import pandas as pd
from collections import defaultdict

from ..pafstats import parse_paf
from ..config import Config, ParamGroup
from ..argparse import ArgParser, Opt 
from ..index import BWA_OPTS
from ..fast5 import Fast5Reader, FAST5_OPTS
from ..sigproc import ProcRead
from .. import DTWd, DTWp, StaticBDTW, BandedDTW, DTW_GLOB, nt

from _uncalled._nt import KmerArray

from .. import PoreModel
from ..dtw import BcFast5Aln, ReadAln, AlnTrack, RefCoord

#TODO make this better
METHODS = {
    "GuidedBDTW" : BandedDTW,
    "StaticBDTW" : StaticBDTW,
    "DTW" : DTWd,
}

class AlignParams(ParamGroup):
    _name = "dtw"

AlignParams._def_params(
    ("method", "GuidedBDTW", str, "DTW method"),
    ("band_width", 50, int, "DTW band width (only applies to BDTW)"),
    ("band_shift", 0.5, float, "DTW band shift coefficent (only applies to BDTW)"),
    #("ref_bounds", None, RefCoord, "Will only output DTW within these reference coordinates if specified"),
    ("mm2_paf", None, str, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string."),
    ("out_path", None, str, "Path to directory where alignments will be stored. If not specified will display interactive dotplot for each read."),
)

OPTS = (Opt("index_prefix", "track"),) + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "dtw", required=True),
    Opt(("-o", "--out-path"), "dtw"),
    Opt("--name", "track"),
    Opt(("-f", "--overwrite"), "track", action="store_true", help="Will overwrite alignment track if one already exists"),
    Opt("--rna", fn="set_r94_rna"),
    Opt(("-R", "--ref-bounds"), "track"),
    Opt("--method", "dtw", choices=METHODS.keys()),
    Opt(("-b", "--band-width"), "dtw"),
    Opt(("-s", "--band-shift"), "dtw"),
    Opt(("-N", "--norm-len"), "normalizer", "len", default=0),
)

def main(conf):
    """Perform DTW alignment guided by basecalled alignments"""
    conf.fast5_reader.load_bc = True
    conf.proc_read.detect_events = True
    conf.export_static()

    fast5s = Fast5Processor(conf=conf)

    track = AlnTrack(conf.dtw.out_path, mode="w", conf=conf)

    mm2s = {p.qr_name : p
         for p in parse_paf(
            conf.dtw.mm2_paf,
            ref_bounds=conf.track.ref_bounds,
            full_overlap=conf.track.full_overlap,
    )}

    t = time.time()
    for read in fast5s:
        db1_time = db2_time = dtw_time = 0
        t = time.time()

        paf = mm2s.get(read.id, None)

        if paf is None:
            continue

        ref_bounds = RefCoord(paf.rf_name, paf.rf_st, paf.rf_en, paf.is_fwd)
        track.init_read_aln(read, ref_bounds)

        db1_time += time.time() - t
        t = time.time()

        if track.read_aln is None:
            sys.stderr.write("Should not happen")
            continue

        dtw = GuidedDTW(track, read, paf, conf)

        if dtw.empty:
            sys.stderr.write("Should not happen?")
            continue

        dtw_time += time.time() - t
        t = time.time()

        track.save_aln()

        db2_time += time.time() - t
        t = time.time()
        
        print("%s\t%f\t%f\t%f" % (read.id, dtw_time, db1_time, db2_time))

    if not track is None:
        t = time.time()
        track.close()
        print("end\t0\t0\t%f" % (time.time()-t))

class GuidedDTW:

    #TODO do more in constructor using prms, not in main
    def __init__(self, track, read, paf, conf=None, **kwargs):
        self.conf = read.conf if conf is None else conf
        self.prms = self.conf.dtw

        self.track = track

        #paf_st, paf_en = self.track.index.ref_to_mref(
        #    paf.rf_name, paf.rf_st, paf.rf_en, paf.is_fwd, self.conf.is_rna)
        #mrefs = pd.RangeIndex(paf_st+nt.K-1, paf_en)

        self.ref_kmers = self.track.load_aln_kmers(store=False)

        self.bcaln = BcFast5Aln(self.track.index, read, paf, self.track.read_aln.mrefs)
        if self.bcaln.empty:
            self.empty = True
            return

        if self.bcaln.errs is not None:
            bcerr = self.bcaln.errs[["mref", "type", "seq"]]#.dropna()
            bcerr.set_index("mref", drop=True, inplace=True)
            bcerr["type"].astype("category", copy=False)
            bcerr = bcerr[~bcerr.index.duplicated(keep="first")]
            self.track.read_aln.set_bcerr(bcerr.sort_index())

        self.read = read

        self.method = self.prms.method
        if not self.method in METHODS:
            sys.stderr.write("Error: unrecongized DTW method \"%s\".\n" % method)
            sys.stderr.write("Must be one of \"%s\".\n" % "\", \"".join(METHODS.keys()))
            sys.exit()

        self.dtw_fn = METHODS[self.method]

        self.model = PoreModel(self.conf.pore_model)

        self.samp_min = self.bcaln.aln['sample'].min()
        self.samp_max = self.bcaln.aln['sample'].max()

        self._calc_dtw()

        self.empty = False

    #TODO refactor inner loop to call function per-block
    def _calc_dtw(self):
        self.mats = list()

        path_qrys = list()
        path_refs = list()

        band_blocks = list()

        bc = self.bcaln.aln#.sort_values("sample").reset_index()
        block_min = int(bc.index[bc['sample'].searchsorted(self.samp_min)])
        block_max = int(bc.index[bc['sample'].searchsorted(self.samp_max)])

        block_starts = np.insert(self.bcaln.ref_gaps, 0, block_min)
        block_ends   = np.append(self.bcaln.ref_gaps, block_max)

        #TODO make this actually do something for spliced RNA
        for st, en in [(block_min, block_max)]:
        #for st, en in zip(block_starts, block_ends):
            samp_st = bc.loc[st,'sample']
            samp_en = bc.loc[en-1,'sample']

            mref_st = st#+nt.K-1#bc.loc[st,"mref"]+nt.K-1
            mref_en = en+1#bc.loc[en-1,"mref"]+1

            read_block = self.read.sample_range(samp_st, samp_en)

            block_signal = read_block['norm_sig'].to_numpy()
            block_kmers = self.ref_kmers.loc[mref_st:mref_en]

            args = self._get_dtw_args(bc, read_block, mref_st, block_kmers)

            dtw = self.dtw_fn(*args)

            path = np.flip(dtw.path)
            #TODO shouldn't need to clip, error in bdtw
            path_qrys.append(read_block.index[np.clip(path['qry'], 0, len(read_block))])
            path_refs.append(mref_st + path['ref'])

            if hasattr(dtw, "ll"):
                band_blocks.append(
                    self._ll_to_df(dtw.ll, read_block, mref_st, len(block_kmers))
                )

        df = pd.DataFrame({'mref': np.concatenate(path_refs)}, 
                               index = np.concatenate(path_qrys),
                               dtype='Int32') \
                  .join(self.read.df) \
                  .drop(columns=['mean', 'stdv', 'mask'], errors='ignore') \
                  .rename(columns={'norm_sig' : 'current'})

        #df['kmer'] = self.ref_kmers.loc[df['mref']].to_numpy()

        self.track.read_aln.set_subevent_aln(df, True)

        if len(band_blocks) == 0:
            self.track.read_aln.bands = None
        elif len(band_blocks) > 1:
            self.track.read_aln.bands = pd.concat(band_blocks)
        else:
            self.track.read_aln.bands = pd.DataFrame(band_blocks[0])

    def _get_dtw_args(self, bc, read_block, mref_start, ref_kmers):
        common = (
            read_block['norm_sig'].to_numpy(), 
            nt.kmer_array(ref_kmers), 
            self.model)
        qry_len = len(read_block)
        ref_len = len(ref_kmers)

        #should maybe move to C++
        if self.method == "GuidedBDTW":
            band_count = qry_len + ref_len
            band_lls = list()

            starts = bc.index[bc['sample'].searchsorted(read_block['start'])]

            q = r = 0
            shift = int(np.round(self.prms.band_shift*self.prms.band_width))
            for i in range(band_count):
                band_lls.append( (int(q+shift), int(r-shift)) )

                tgt = starts[q] if q < len(starts) else starts[-1]
                if r <= tgt - mref_start:
                    r += 1
                else:
                    q += 1

            return common + (self.prms.band_width, band_lls)

        elif self.method == "StaticBDTW":
            return common + (self.prms.band_width, self.prms.band_shift)

        else:
            return common + (DTW_GLOB,)

    #TODO store in ReadAln metadata
    def _ll_to_df(self, ll, read_block, ref_st, ref_len):
        block_qry_st = np.clip(ll['qry'],                 0, len(read_block)-1)
        block_qry_en = np.clip(ll['qry']-self.prms.band_width, 0, len(read_block)-1)
        block_ref_st = np.clip(ll['ref'],                 0, ref_len-1)
        block_ref_en = np.clip(ll['ref']+self.prms.band_width, 0, ref_len-1)
                                              
        band_samps = read_block['start'].to_numpy()
        band_ref_st = np.zeros(len(band_samps))
        band_ref_en = band_ref_st + ref_len

        band_ref_st[block_qry_st] = block_ref_st
        band_ref_en[block_qry_en] = block_ref_en+1

        return pd.DataFrame({
            'samp': band_samps, 
            'ref_st': ref_st + band_ref_st, 
            'ref_en': ref_st + band_ref_en
        })


class Fast5Processor(Fast5Reader):
    def __next__(self):
        return ProcRead(Fast5Reader.__next__(self), conf=self.conf)
    
    def __getitem__(self, read_id):
        return ProcRead(Fast5Reader.__getitem__(self, read_id), conf=self.conf)
