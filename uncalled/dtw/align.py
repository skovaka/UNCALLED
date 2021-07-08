import sys
import time
import re
import numpy as np
import pandas as pd

from ..pafstats import parse_paf
from ..config import Config, ArgParser, ParamGroup, Opt
from ..index import BWA_OPTS
from ..fast5 import Fast5Reader, FAST5_OPTS
from ..sigproc import ProcRead
from .. import BwaIndex, DTWd, DTWp, StaticBDTW, BandedDTW, DTW_GLOB, nt

from .dotplot import Dotplot
from . import PoreModel, BcFast5Aln, ReadAln, Track, ref_coords

#TODO make this better
METHODS = {
    "GuidedBDTW" : BandedDTW,
    "StaticBDTW" : StaticBDTW,
    "DTW" : DTWd,
}

class AlignParams(ParamGroup):
    _name = "align"

AlignParams._def_params(
    ("method", "GuidedBDTW", str, "DTW method"),
    ("band_width", 50, int, "DTW band width (only applies to BDTW)"),
    ("band_shift", 0.5, float, "DTW band shift coefficent (only applies to BDTW)"),
    ("ref_bounds", None, tuple, "Will only output DTW within these reference coordinates if specified"),
    ("mm2_paf", None, str, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string."),
    ("out_path", None, str, "Path to directory where alignments will be stored. If not specified will display interactive dotplot for each read."),
)

OPTS = BWA_OPTS + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "align", required=True),
    Opt(("-o", "--out-path"), "align"),
    Opt(("-f", "--overwrite"), "track", action="store_true", help="Will overwrite alignment track if one already exists"),
    Opt("--rna", fn="set_r94_rna"),
    Opt(("-R", "--ref-bounds"), "align", type=ref_coords),
    Opt("--method", "align", choices=METHODS.keys()),
    Opt(("-b", "--band-width"), "align"),
    Opt(("-s", "--band-shift"), "align"),
    Opt(("-N", "--norm-len"), "normalizer", "len", default=0),
)

def main(conf):
    """Performs DTW alignment and either outputs as alignment Track or displays dotplots"""
    conf.fast5_reader.load_bc = True
    conf.proc_read.detect_events = True
    conf.export_static()

    mm2s = dict()
    if conf.align.mm2_paf is not None:
        for p in parse_paf(conf.align.mm2_paf):
            old = mm2s.get(p.qr_name, None)
            if old is None or old.aln_len < p.aln_len:
                mm2s[p.qr_name] = p

    idx = BwaIndex(conf.mapper.bwa_prefix, True)

    fast5s = Fast5Processor(conf=conf)

    if conf.align.out_path is not None:
        track = Track(conf.align.out_path, mode="w", conf=conf)
    else:
        track = None

    t = time.time()

    for read in fast5s:

        if not read.id in mm2s:
            continue

        dtw = GuidedDTW(idx, read, mm2s[read.id], conf)

        if dtw.empty:
            continue

        if track is None:
            dplt = Dotplot(idx, read, conf=conf)
            dplt.add_aln(dtw.bcaln, False)
            dplt.add_aln(dtw.aln, True)
            dplt.show()
        else:
            track.save_aln(dtw.aln, read.f5.filename)
        
        print(read.id)

class GuidedDTW:

    #TODO do more in constructor using prms, not in main
    def __init__(self, index, read, paf, conf=None, **kwargs):
        self.conf = read.conf if conf is None else conf
        self.prms = self.conf.align

        t = time.time()

        self.bcaln = BcFast5Aln(index, read, paf, self.conf.align.ref_bounds)
        if self.bcaln.empty:
            self.empty = True
            return

        self.aln = ReadAln(index, paf, is_rna=self.conf.is_rna)

        self.read = read
        self.idx = index

        self.method = self.prms.method
        if not self.method in METHODS:
            sys.stderr.write("Error: unrecongized DTW method \"%s\".\n" % method)
            sys.stderr.write("Must be one of \"%s\".\n" % "\", \"".join(METHODS.keys()))
            sys.exit()

        self.dtw_fn = METHODS[self.method]

        self.model = PoreModel(self.conf.pore_model)

        self.ref_name = self.bcaln.ref_name
        self.ref_min = self.bcaln.ref_start
        self.ref_max = self.bcaln.ref_end

        self.samp_min = self.bcaln.df['sample'].min()
        self.samp_max = self.bcaln.df['sample'].max()

        self.ref_kmers = self.aln.get_index_kmers(self.idx)
        self._calc_dtw()

        self.empty = False

    #TODO refactor inner loop to call function per-block
    def _calc_dtw(self):
        self.mats = list()

        path_qrys = list()
        path_refs = list()

        band_blocks = list()

        block_min = self.bcaln.df['sample'].searchsorted(self.samp_min)
        block_max = self.bcaln.df['sample'].searchsorted(self.samp_max)


        block_starts = np.insert(self.bcaln.ref_gaps, 0, block_min)
        block_ends   = np.append(self.bcaln.ref_gaps, block_max)

        #TODO make this actually do something for spliced RNA
        for st, en in [(block_min, block_max)]:
        #for st, en in zip(block_starts, block_ends):
            samp_st = self.bcaln.df.loc[st,'sample']
            samp_en = self.bcaln.df.loc[en-1,'sample']

            refmir_st = self.bcaln.df.loc[st,"refmir"]
            refmir_en = self.bcaln.df.loc[en-1,"refmir"]

            #print(samp_st, samp_en, refmir_st, refmir_en, self.aln.refmir_start, self.idx.refmir_to_ref(refmir_en))

            read_block = self.read.sample_range(samp_st, samp_en)

            block_signal = read_block['norm_sig'].to_numpy()
            block_kmers = self.ref_kmers[refmir_st-self.aln.refmir_start:refmir_en-self.aln.refmir_start]

            args = self._get_dtw_args(read_block, refmir_st, block_kmers)

            dtw = self.dtw_fn(*args)

            #TODO flip in traceback
            path = np.flip(dtw.path)
            path_qrys.append(read_block.index[path['qry']])
            path_refs.append(refmir_st + path['ref'])

            if hasattr(dtw, "ll"):
                band_blocks.append(
                    self._ll_to_df(dtw.ll, read_block, refmir_st, len(block_kmers))
                )

        df = pd.DataFrame({'refmir': np.concatenate(path_refs)}, 
                               index = np.concatenate(path_qrys),
                               dtype='Int32') \
                  .join(self.read.df) \
                  .drop(columns=['mean', 'stdv', 'mask'], errors='ignore') \
                  .rename(columns={'norm_sig' : 'mean'})

        y_min = self.aln.refmir_start
        df['kmer'] = self.ref_kmers[df['refmir'].astype(int)-y_min]

        self.aln.set_subevent_aln(df, True)

        #self.dtw['refmir'] += self.bcaln.y_min

        if len(band_blocks) == 0:
            self.aln.bands = None
        elif len(band_blocks) > 1:
            self.aln.bands = pd.concat(band_blocks)
        else:
            self.aln.bands = pd.DataFrame(band_blocks[0])

    def _get_dtw_args(self, read_block, ref_start, ref_kmers):
        common = (read_block['norm_sig'].to_numpy(), ref_kmers, self.model)
        qry_len = len(read_block)
        ref_len = len(ref_kmers)

        #should maybe move to C++
        if self.method == "GuidedBDTW":
            band_count = qry_len + ref_len
            band_lls = list()

            starts = self.bcaln.df['sample'].searchsorted(read_block['start'])

            q = r = 0
            shift = int(np.round(self.prms.band_shift*self.prms.band_width))
            for i in range(band_count):
                band_lls.append( (int(q+shift), int(r-shift)) )

                tgt = starts[q] if q < len(starts) else starts[-1]
                if r <= self.bcaln.df.loc[tgt,'refmir'] - ref_start:
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
