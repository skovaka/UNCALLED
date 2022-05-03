import sys
import time
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import progressbar as progbar

from sklearn.linear_model import TheilSenRegressor
from ..pafstats import parse_paf
from ..config import Config, ParamGroup
from ..argparse import ArgParser, Opt 
from ..index import BWA_OPTS, str_to_coord
from ..fast5 import Fast5Reader, FAST5_OPTS

from _uncalled import (
    GlobalDTWK5, StaticBDTWK5, BandedDTWK5, 
    GlobalDTWK10, StaticBDTWK10, BandedDTWK10, 
    DTW_PRMS_EVT_GLOB, DtwParams
)

from ..signal_processor import SignalProcessor

from .. import PoreModel
from . import Bcaln, Tracks

#TODO make this better
METHODS = {
    "guided" : {5: BandedDTWK5, 10: BandedDTWK10},
    "static" : {5: StaticBDTWK5, 10: StaticBDTWK10},
    "global" : {5: GlobalDTWK5, 10: GlobalDTWK10}
}

#class DtwParams(ParamGroup):
#    _name = "dtw"
#
#DtwParams._def_params(
#    ("method", "guided", str, "DTW method"),
#    ("band_width", 50, int, "DTW band width (only applies to BDTW)"),
#    ("band_shift", 0.5, float, "DTW band shift coefficent (only applies to BDTW)"),
#    ("mm2_paf", None, str, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string."),
##    ("mask_skips", False, bool, "Represent skips as missing data"),
#)

OPTS = (Opt("index_prefix", "tracks"),) + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "dtw", required=True),
    Opt(("-o", "--db-out"), "tracks.io"),
    Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("-a", "--append"), "tracks.io", action="store_true"),
    Opt(("-p", "--pore-model"), "pore_model", "name"),
    Opt("--full-overlap", "tracks", action="store_true"),
    #Opt(("-S", "--mask-skips"), "dtw", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    #Opt("--method", "dtw", choices=METHODS.keys()),
    Opt(("-i", "--iterations"), "dtw"),
    Opt(("-c", "--cost-fn"), "dtw", choices=["abs_diff","z_score","norm_pdf"]),
    Opt("--skip-cost", "dtw"),
    Opt("--stay-cost", "dtw"),
    Opt("--move-cost", "dtw"),
    Opt(("-b", "--band-width"), "dtw"),
    Opt(("-s", "--band-shift"), "dtw"),
    Opt(("-N", "--norm-len"), "normalizer", "len", default=0),
)

def main(conf):
    """Perform DTW alignment guided by basecalled alignments"""
    conf.fast5_reader.load_bc = True
    conf.proc_read.detect_events = True
    conf.export_static()

    tracks = Tracks(conf=conf)

    clip_coords = tracks.coords

    fast5s = Fast5Reader(conf=conf)

    read_filter = fast5s.get_read_filter()

    pafs = parse_paf(
        conf.dtw.mm2_paf, 
        ref_bounds=conf.tracks.ref_bounds, 
        read_filter=read_filter,
        full_overlap=conf.tracks.full_overlap)

    mm2s = defaultdict(list)
    for paf in pafs:
        mm2s[paf.qr_name].append(paf)

    model = PoreModel(conf.pore_model)
    sigproc = SignalProcessor(model, conf)

    #pbar = progbar.ProgressBar(
    #        widgets=[progbar.Percentage(), progbar.Bar(), progbar.ETA()], 
    #        maxval=len(mm2s)).start()

    n_reads = 0

    for read in fast5s:
        aligned = False
        for paf in mm2s[read.id]:
            t0 = time.time()
            dtw = GuidedDTW(tracks, sigproc, read, paf, conf)

            sys.stderr.write(f"{read.id}\n")

            if dtw.df is None:
                sys.stderr.write("# dtw failed\n")
                continue
            aligned = True

        if aligned:
            #pbar.update(n_reads)
            n_reads += 1

    tracks.close()

    #pbar.finish()

class GuidedDTW:

    #TODO do more in constructor using prms, not in main
    #def __init__(self, track, read, paf, conf=None, **kwargs):
    def __init__(self, tracks, sigproc, read, paf, conf=None, **kwargs):
        self.conf = read.conf if conf is None else conf
        self.prms = self.conf.dtw

        bcaln = Bcaln(conf, tracks.index, read, paf, tracks.coords)
        if bcaln.empty:
            self.df = None
            return


        aln_id, self.coords = tracks.write_alignment(read.id, read.filename, bcaln.coords, {"bcaln" : bcaln.df})
        #TODO return coords?

        kmers = self.coords.kmers.sort_index()

        self.bcaln = bcaln.df[bcaln.df.index.isin(kmers.index)].sort_index()[["start"]].dropna()

        self.ref_gaps = list(sorted(bcaln.ref_gaps))

        k = kmers.index[0]
        kmer_blocks = list()
        for start,end in self.ref_gaps:
            kmer_blocks.append(kmers.loc[k:start])
            k = end
        kmer_blocks.append(kmers.loc[k:])

        self.ref_kmers = pd.concat(kmer_blocks)
        #self.ref_kmers = kmers
        

        self.model = PoreModel(self.conf.pore_model)

        self.method = self.prms.band_mode
        if not self.method in METHODS:
            opts = "\", \"".join(METHODS.keys())
            raise ValueError(f"Error: unrecongized DTW method \"{method}. Must be one of \"{opts}\"")

        if not self.model.K in METHODS[self.method]:
            raise ValueError(f"Invalid DTW k-mer length: {self.model.K}")

        self.dtw_fn = METHODS[self.method][self.model.K]

        self.samp_min = self.bcaln["start"].min()
        self.samp_max = self.bcaln["start"].max()

        signal = sigproc.process(read)

        df = self._calc_dtw(signal)

        for i in range(self.prms.iterations-1):
            reg = self.renormalize(signal, df)
            signal.rescale(reg.coef_, reg.intercept_)
            df = self._calc_dtw(signal)

        self.df = df.set_index("mref")

        tracks.write_dtw_events(self.df, aln_id=aln_id)

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

        band_blocks = list()

        #bc = self.bcaln.df.sort_index()
        block_min = self.bcaln.index[0]
        block_max = self.bcaln.index[-1]

        #TODO for spliced RNA, must find gaps in mref index
        #block_starts = np.insert(self.bcaln.ref_gaps, 0, block_min)
        #block_ends   = np.append(self.bcaln.ref_gaps, block_max)
        #for st, en in zip(block_starts, block_ends):

        for st, en in [(block_min, block_max)]:
            mref_st = st
            mref_en = en+1
            block_kmers = self.ref_kmers.loc[mref_st:mref_en]

            samp_st = self.bcaln.loc[st,"start"]
            samp_en = self.bcaln.loc[en,"start"]

            read_block = signal.sample_range(samp_st, samp_en)

            block_signal = read_block['mean'].to_numpy()
            
            args = self._get_dtw_args(read_block, mref_st, block_kmers)

            dtw = self.dtw_fn(*args)
            
            path = np.flip(dtw.path)
            #TODO shouldn't need to clip, error in bdtw
            path_qrys.append(read_block.index[np.clip(path['qry'], 0, len(read_block))])
            path_refs.append(block_kmers.index[path['ref']])

            if hasattr(dtw, "ll"):
                band_blocks.append(
                    self._ll_to_df(dtw.ll, read_block, mref_st, len(block_kmers))
                )

        return pd.DataFrame({'mref': np.concatenate(path_refs)}, 
                               index = np.concatenate(path_qrys),
                               dtype='Int32') \
                  .join(signal.to_df()) \
                  .drop(columns=['mask'], errors='ignore') \
                  .rename(columns={'mean' : 'current', 'stdv' : 'current_stdv'})

    def _get_dtw_args(self, read_block, mref_start, ref_kmers):
        common = (
            self.prms,
            read_block['mean'].to_numpy(), 
            self.model.kmer_array(ref_kmers),  #TODO probably dont need convert
            self.model.instance)
        qry_len = len(read_block)
        ref_len = len(ref_kmers)

        #should maybe move to C++
        if self.method == "guided":
            band_count = qry_len + ref_len
            band_lls = list()

            starts = self.bcaln.index[self.bcaln["start"].searchsorted(read_block['start'])]

            q = r = 0
            shift = int(np.round(self.prms.band_shift*self.prms.band_width))
            for i in range(band_count):
                band_lls.append( (int(q+shift), int(r-shift)) )
                #print(band_lls[-1])

                tgt = starts[q] if q < len(starts) else starts[-1]
                if r < len(self.ref_kmers) and self.ref_kmers.index[r] <= tgt:
                    r += 1
                else:
                    q += 1

            return common + (band_lls, )

        elif self.method == "static":
            return common

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

class Fast5Processor(Fast5Reader):
    def __next__(self):
        return ProcRead(Fast5Reader.__next__(self), conf=self.conf)
    
    def __getitem__(self, read_id):
        return ProcRead(Fast5Reader.__getitem__(self, read_id), conf=self.conf)
