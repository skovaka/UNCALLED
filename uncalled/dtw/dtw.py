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
from . import Bcaln, AlnTrack, RefCoord, TrackIO

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

OPTS = (Opt("index_prefix", "track_io"),) + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "dtw", required=True),
    Opt(("-o", "--out-path"), "dtw"),
    Opt("--name", "track"),
    Opt(("-f", "--overwrite"), "track", action="store_true", help="Will overwrite alignment track if one already exists"),
    Opt("--rna", fn="set_r94_rna"),
    Opt(("-R", "--ref-bounds"), "track_io"),
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

    #track = AlnTrack(conf.dtw.out_path, mode="w", conf=conf, load_mat=False)
    track_io = TrackIO(None, conf.dtw.out_path, conf=conf)

    mm2s = {p.qr_name : p
         for p in parse_paf(
            conf.dtw.mm2_paf,
            ref_bounds=conf.track_io.ref_bounds,
            full_overlap=conf.track_io.full_overlap,
    )}

    for read in fast5s:

        paf = mm2s.get(read.id, None)

        if paf is None:
            continue

        print(read.id)
        dtw = GuidedDTW(track_io, read, paf, conf)

        if dtw.df is None:
            sys.stderr.write("dtw failed\n")
            continue

    track_io.close()

class GuidedDTW:

    #TODO do more in constructor using prms, not in main
    #def __init__(self, track, read, paf, conf=None, **kwargs):
    def __init__(self, track_io, read, paf, conf=None, **kwargs):
        self.conf = read.conf if conf is None else conf
        self.prms = self.conf.dtw

        #self.track = track

        bcaln = Bcaln(track_io.index, read, paf, track_io.coords)
        if bcaln.empty:
            self.df = None
            return


        #TODO init_alignment(read_id, fast5, group, layers)
        self.track = track_io.init_alignment(read.id, read.f5.filename, "bcaln", bcaln.df)

        self.bcaln = bcaln.df.sort_index()

        #self.ref_kmers = self.track.load_aln_kmers().sort_index()
        #print(self.ref_kmers)
        self.ref_kmers = self.track.coords.load_kmers(track_io.index).sort_index()

        self.read = read

        self.method = self.prms.method
        if not self.method in METHODS:
            sys.stderr.write("Error: unrecongized DTW method \"%s\".\n" % method)
            sys.stderr.write("Must be one of \"%s\".\n" % "\", \"".join(METHODS.keys()))
            sys.exit()

        self.dtw_fn = METHODS[self.method]

        self.model = PoreModel(self.conf.pore_model)

        self.samp_min = self.bcaln['sample'].min()
        self.samp_max = self.bcaln['sample'].max()

        self._calc_dtw()

        self.empty = False

    #TODO refactor inner loop to call function per-block
    def _calc_dtw(self):
        self.mats = list()

        path_qrys = list()
        path_refs = list()

        band_blocks = list()

        #bc = self.bcaln.df.sort_index()
        block_min = self.bcaln.index[0]
        block_max = self.bcaln.index[-1]

        #TODO for spliced RNA, must find gaps in mref index
        #block_starts = np.insert(self.bcaln.ref_gaps, 0, block_min)
        #block_ends   = np.append(self.bcaln.ref_gaps, block_max)
        #for st, en in zip(block_starts, block_ends):

        for st, en in [(block_min, block_max)]:
            samp_st = self.bcaln.loc[st,'sample']
            samp_en = self.bcaln.loc[en,'sample']

            mref_st = st
            mref_en = en+1

            read_block = self.read.sample_range(samp_st, samp_en)

            block_signal = read_block['norm_sig'].to_numpy()
            block_kmers = self.ref_kmers.loc[mref_st:mref_en]

            args = self._get_dtw_args(read_block, mref_st, block_kmers)

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

        #def self, aln_id, group, layers):
        self.df = collapse_events(df, True)
        self.track.add_layer_group("dtw", self.df)

        #if len(band_blocks) == 0:
        #    self.track.read_aln.bands = None
        #elif len(band_blocks) > 1:
        #    self.track.read_aln.bands = pd.concat(band_blocks)
        #else:
        #    self.track.read_aln.bands = pd.DataFrame(band_blocks[0])

    def _get_dtw_args(self, read_block, mref_start, ref_kmers):
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

            starts = self.bcaln.index[self.bcaln['sample'].searchsorted(read_block['start'])]

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

def collapse_events(dtw, kmer_str=False, start_col="start", length_col="length", mean_col="current", kmer_col="kmer"):

    dtw["cuml_mean"] = dtw[length_col] * dtw[mean_col]

    grp = dtw.groupby("mref")

    if kmer_col in dtw:
        if kmer_str:
            kmers = [nt.kmer_rev(nt.str_to_kmer(k,0)) for k in grp[kmer_col].first()]
        else:
            kmers = grp[kmer_col].first()
    else:
        kmers = None

    mrefs = grp["mref"].first()

    lengths = grp[length_col].sum()

    dtw = pd.DataFrame({
        "mref"    : mrefs.astype("int64"),
        "start"  : grp[start_col].min().astype("uint32"),
        "length" : lengths.astype("uint32"),
        "current"   : grp["cuml_mean"].sum() / lengths
    })

    if kmers is not None:
        dtw["kmer"] = kmers.astype("uint16")

    return dtw.set_index("mref").sort_index()

class Fast5Processor(Fast5Reader):
    def __next__(self):
        return ProcRead(Fast5Reader.__next__(self), conf=self.conf)
    
    def __getitem__(self, read_id):
        return ProcRead(Fast5Reader.__getitem__(self, read_id), conf=self.conf)
