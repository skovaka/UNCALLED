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
from _uncalled import PORE_MODELS, BwaIndex, DTWd, DTWp, StaticBDTW, BandedDTW, DTW_GLOB, nt

from .dotplot import Dotplot
from . import BcFast5Aln, ReadAln, Track, ref_coords

#TODO make this better
METHODS = {
    "DTWd" : DTWd,
    "DTWp" : DTWp,
    "StaticBDTW" : StaticBDTW,
    "GuidedBDTW" : BandedDTW
}

class AlignParams(ParamGroup):
    _name = "align"

AlignParams._def_params(
    ("method", "DTWd", str, "DTW method"),
    ("band_width", 100, int, "DTW band width (only applies to BDTW)"),
    ("band_shift", 0.5, float, "DTW band shift coefficent (only applies to BDTW)"),
    ("ref_bounds", None, tuple, "Will only output DTW within these reference coordinates if specified"),
    ("mm2_paf", None, str, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string."),
    ("out_path", None, str, "Path to directory where alignments will be stored. If not specified will display interactive dotplot for each read."),
)

OPTS = BWA_OPTS + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "align", required=True),
    Opt(("-c", "--start-chunk"), "read_buffer"),
    Opt(("-C", "--max-chunks"), "read_buffer"),
    Opt("--rna", fn="set_r94_rna"),
    Opt("--method", "align", choices=METHODS.keys()),
    Opt(("-b", "--band-width"), "align"),
    Opt(("-s", "--band-shift"), "align"),
    Opt(("-N", "--norm-len"), "normalizer", "len", default=0),
    Opt(("-R", "--ref-bounds"), "align", type=ref_coords),
    Opt(("-f", "--force-overwrite"), action="store_true"),
    Opt(("-o", "--out-path"), "align"),
)

def main(conf):
    """Performs DTW alignment and either outputs as alignment Track or displays dotplots"""

    #if conf.rna:
    #    conf.set_r94_rna()

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
        track = Track(conf.align.out_path, "w", conf=conf, overwrite=conf.force_overwrite)
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
    def __init__(self, index, read, paf, conf=None, dtw_events=None, **kwargs):
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

        model_name = self.conf.mapper.pore_model

        #TODO clean this up
        if model_name.endswith("_compl"):
            model_name = model_name[:-5]+"templ"

        self.model = PORE_MODELS[model_name]

        self.ref_name = self.bcaln.ref_name
        self.ref_min = self.bcaln.ref_start
        self.ref_max = self.bcaln.ref_end

        self.samp_min = self.bcaln.df['sample'].min()
        self.samp_max = self.bcaln.df['sample'].max()

        self.ref_kmers = self.aln.get_index_kmers(self.idx)

        if dtw_events is None:
            self.calc_dtw()
            #self.calc_events()
        else:
            self.load_dtw_events(dtw_events)

        self.empty = False

    #TODO generate AlignedRead
    def calc_events(self):

        if self.read.has_events:
            self.aln.df['sum'] = self.aln.df['signal'] * self.aln.df['length']

        grp = self.aln.df.groupby("miref")
        sigs = grp['signal']

        ref_coords = np.abs(self.bcaln.y_min + grp['miref'].first())

        self.events = pd.DataFrame({
            "miref"   : ref_coords,
            "start"  : grp['sample'].min(),
            "kmer" : grp['kmer'].first(),
        })#.reset_index(drop=True)

        if self.read.has_events:
            self.events['length'] = grp['length'].sum()
            self.events['mean'] = grp['sum'].sum() / self.events['length']
            self.aln.df.drop(columns=['sum'])
        else:
            self.events['length'] = grp['signal'].count()
            self.events['mean'] = grp['signal'].mean() 
            self.events['stdv'] = grp['signal'].std().fillna(0)

    #TODO do with mirror coords in BWA index
    def load_kmers(self):

        shift = nt.K - 1

        pac_st, pac_en = self.idx.mirror_ref_coords(self.ref_name, self.ref_min, self.ref_max, self.bcaln.is_fwd, not self.bcaln.seq_fwd)
        pac_st -= shift

        kmers2 = self.idx.get_kmers(pac_st, pac_en, not self.bcaln.seq_fwd)
        kmers3 = self.aln.get_index_kmers(self.idx)

        st = self.ref_min - (shift if not self.bcaln.flip_ref else 0)
        en = self.ref_max + (shift if self.bcaln.flip_ref else 0)

        if st < 0:
            pad = -st
            st = 0
        else:
            pad = 0

        kmers = np.array(self.idx.get_kmers(
            self.bcaln.rf_name, st, en
        ))


        if self.bcaln.flip_ref:
            kmers = np.flip(nt.kmer_rev(kmers))

        if not self.bcaln.is_fwd:
            kmers = nt.kmer_comp(kmers)

    def get_dtw_args(self, read_block, ref_start, ref_kmers):
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
                if r <= self.bcaln.df.loc[tgt,'miref'] - ref_start:
                    r += 1
                else:
                    q += 1

            return common + (self.prms.band_width, band_lls)

        elif self.method == "StaticBDTW":
            return common + (self.prms.band_width, self.prms.band_shift)

        else:
            return common + (DTW_GLOB,)

    #TODO store in ReadAln metadata
    def ll_to_df(self, ll, read_block, ref_st, ref_len):
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


    #TODO refactor inner loop to call function per-block
    def calc_dtw(self):
        self.mats = list()

        path_qrys = list()
        path_refs = list()

        band_blocks = list()

        block_min = self.bcaln.df['sample'].searchsorted(self.samp_min)
        block_max = self.bcaln.df['sample'].searchsorted(self.samp_max)

        y_min = self.aln.miref_start

        block_starts = np.insert(self.bcaln.ref_gaps, 0, block_min)
        block_ends   = np.append(self.bcaln.ref_gaps, block_max)

        #TODO make this actually do something for spliced RNA
        for st, en in [(block_min, block_max)]:
        #for st, en in zip(block_starts, block_ends):
            samp_st = self.bcaln.df.loc[st,'sample']
            samp_en = self.bcaln.df.loc[en-1,'sample']

            miref_st = self.bcaln.df.loc[st,"miref"]
            miref_en = self.bcaln.df.loc[en-1,"miref"]

            read_block = self.read.sample_range(samp_st, samp_en)

            block_signal = read_block['norm_sig'].to_numpy()
            block_kmers = self.ref_kmers[miref_st-self.aln.miref_start:miref_en-self.aln.miref_start]

            args = self.get_dtw_args(read_block, miref_st, block_kmers)

            dtw = self.dtw_fn(*args)

            #TODO flip in traceback
            path = np.flip(dtw.path)
            print(list(dtw.path))
            path_qrys.append(read_block.index[path['qry']])
            path_refs.append(miref_st + path['ref'])

            if hasattr(dtw, "ll"):
                band_blocks.append(
                    self.ll_to_df(dtw.ll, read_block, miref_st, len(block_kmers))
                )

        df = pd.DataFrame({'miref': np.concatenate(path_refs)}, 
                               index = np.concatenate(path_qrys),
                               dtype='Int32') \
                  .join(self.read.df) \
                  .drop(columns=['mean', 'stdv', 'mask'], errors='ignore') \
                  .rename(columns={'norm_sig' : 'mean'})
        df['kmer'] = self.ref_kmers[df['miref'].astype(int)-y_min]

        self.aln.set_subevent_aln(df, True)

        #self.dtw['miref'] += self.bcaln.y_min

        if len(band_blocks) == 0:
            self.aln.bands = None
        elif len(band_blocks) > 1:
            self.aln.bands = pd.concat(band_blocks)
        else:
            self.aln.bands = pd.DataFrame(band_blocks[0])

    #TODO move to AlignedRead
    def load_dtw_events(self, event_file):
        self.events = pd.read_pickle(event_file).reset_index()

        block_min = self.bcaln.df['sample'].searchsorted(self.samp_min)
        y_min1 = self.bcaln.df['miref'][block_min]

        y_min = self.events['miref'].min()
        y_max = self.events['miref'].max()

        self.events = self.events.loc[(self.events['start'] >= self.samp_min) & (self.events['start'] <= self.samp_max)]#.reset_index(drop=True)

        y_min2 = self.events['miref'].min()

        if self.bcaln.flip_ref:
            self.events['idx'] = -self.events['miref'] + y_max
        else:
            self.events['idx'] = self.events['miref'] - y_min

        self.events.set_index('idx', inplace=True)
        self.events.sort_index(inplace=True)

        self.aln.df = self.events.drop(columns=["miref"]) \
                              .reset_index() \
                              .rename(columns={'idx' : 'miref', 'start' : 'sample', 'mean' : 'signal'})
        self.aln.df.reset_index()
        self.bands = None

    #TODO move to dotplot
    def plot_dotplot(self, ax):
        if self.bands is not None:
            ax.fill_between(self.bands['samp'], self.bands['ref_st']-1, self.bands['ref_en'], zorder=1, color='#ccffdd', linewidth=1, edgecolor='black', alpha=0.5)

        return ax.step(self.aln.df['start'], self.aln.df['miref'],where="post",color="purple", zorder=3, linewidth=3)
    
    def plot_signal(self, ax_sig):
        samp_min, samp_max = self.aln.get_samp_bounds()

        samps = np.arange(samp_min, samp_max)
        raw_norm = self.read.get_norm_signal(samp_min, samp_max)

        ymin = np.min(raw_norm[raw_norm>0])
        ymax = np.max(raw_norm[raw_norm>0])
        bases = nt.kmer_base(self.aln.df['kmer'], 2)

        samp_bases = np.zeros(len(samps), int)
        for i in range(len(self.aln.df)):
            st = int(self.aln.df.iloc[i]['start'] - samp_min)
            en = int(st + self.aln.df.iloc[i]['length'])
            samp_bases[st:en] = bases[i]

        BASE_COLORS = [
            "#80ff80",
            "#8080ff",
            "#ffbd00",
            "#ff8080",
        ]
        for base, color in enumerate(BASE_COLORS):
            ax_sig.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color, interpolate=True)

        ax_sig.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, alpha=0.75, c="#777777")

        ax_sig.step(self.aln.df['start'], self.model.get_mean(self.aln.df['kmer']), color='white', linewidth=2, where="post")

        ax_sig.vlines(self.aln.df['start'], ymin, ymax, linewidth=2, color="white")

        evts = (self.read.df['start'] >= samp_min) & (self.read.df['start'] < samp_max) & (self.read.df['norm_sig'] > 0)

        if self.read.has_events:
            ax_sig.step(self.read.df['start'][evts], self.read.df['norm_sig'][evts], where='post', color='black', linewidth=3)
        else:
            ax_sig.scatter(self.read.df['start'][evts], self.read.df['norm_sig'][evts], s=5, alpha=0.75, c="#777777") 

    #TODO move to dotplot
    def plot_dtw_events(self, ax_sig, ax_padiff):
        c = 'purple'

        self.plot_signal(ax_sig)

        model_means = self.model.get_mean(self.aln.df['kmer'])

        pa_diffs = np.abs(self.aln.df['mean'] - self.model.get_mean(self.aln.df['kmer']))

        ax_padiff.step(pa_diffs, self.aln.df['miref'], color=c, where="post")

#TODO move to ReadAln
def save(dtw, track):
    events_out = dtw.aln.df.sort_index()

    track.add_read(dtw.read.id, dtw.read.f5.filename, events_out)

class Fast5Processor(Fast5Reader):
    def __next__(self):
        return ProcRead(Fast5Reader.__next__(self), conf=self.conf)
    
    def __getitem__(self, read_id):
        return ProcRead(Fast5Reader.__getitem__(self, read_id), conf=self.conf)
