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

from .dotplot import dotplot
from . import Track, ref_coords

#from . import T

#TODO make this better
METHODS = {
    "DTWd" : DTWd,
    "DTWp" : DTWp,
    "StaticBDTW" : StaticBDTW,
    "GuidedBDTW" : BandedDTW
}

class AlignParams(ParamGroup): pass
AlignParams._def_params(
    ("method", "DTWd", str, "DTW method"),
    ("band_width", 100, int, "DTW band width (only applies to BDTW)"),
    ("band_shift", 0.5, float, "DTW band shift coefficent (only applies to BDTW)"),
    ("ref_bounds", None, tuple, "Will only output DTW within these reference coordinates if specified"),
    ("mm2_paf", None, str, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string."),
    ("out_path", None, str, "Path to directory where alignments will be stored. If not specified will display interactive dotplot for each read."),
)
Config._EXTRA_GROUPS["align"] = AlignParams #TODO put in ParamGroup con

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

    print(type(conf.align.mm2_paf), conf.align.mm2_paf)

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

    for read in fast5s:
        bcaln = BcFast5Aln(read, mm2s.get(read.id, None))
        if bcaln.empty:
            continue

        t = time.time()

        #TODO move to dotplot main
        #if track is not None:
        #    aln_file = track.aln_fname(read.id)
        #else:
        #    aln_file = None

        dtw = GuidedDTW(idx, read, bcaln, conf)

        if dtw.empty:
            continue

        print(read.id)
        if track is None:
            dotplot(dtw)
        else:
            save(dtw, track)

            #TODO add progressbar

class GuidedDTW:

    #TODO do more in constructor using prms, not in main
    def __init__(self, index, read, bcaln, conf=None, dtw_events=None, **kwargs):
        self.conf = read.conf if conf is None else conf
        self.prms = self.conf.align

        #self.prms.set(**kwargs)

        t = time.time()

        self.read = read
        self.bcaln = bcaln
        self.idx = index

        self.method = self.prms.method
        if not self.method in METHODS:
            sys.stderr.write("Error: unrecongized DTW method \"%s\".\n" % method)
            sys.stderr.write("Must be one of \"%s\".\n" % "\", \"".join(METHODS.keys()))
            sys.exit()

        self.dtw_fn = METHODS[self.method]

        model_name = self.conf.mapper.pore_model

        #TODO probably need to rethink fwd/rev compl, but either way clean this up
        if model_name.endswith("_compl"):
            model_name = model_name[:-5]+"templ"

        self.model = PORE_MODELS[model_name]

        self.ref_min = self.bcaln.rf_st
        self.ref_max = self.bcaln.rf_en

        #self.ref_min, self.ref_max = sorted(
        #        np.abs([self.bcaln.y_min, self.bcaln.y_max])
        #)

        if self.prms.ref_bounds is None:
            self.samp_min = self.bcaln.df['sample'].min()
            self.samp_max = self.bcaln.df['sample'].max()
        else:
            self.ref_min = max(self.ref_min, self.prms.ref_bounds[1])
            self.ref_max = min(self.ref_max, self.prms.ref_bounds[2])

            if self.ref_min >= self.ref_max:
                self.empty = True
                return


            self.samp_min = int(bcaln.ref_to_samp(self.ref_min))
            self.samp_max = int(bcaln.ref_to_samp(self.ref_max))
            if bcaln.flip_ref:
                self.samp_min, self.samp_max = (self.samp_max, self.samp_min)


        self.load_kmers()

        if dtw_events is None:
            self.calc_dtw()
            self.calc_events()
        else:
            self.load_dtw_events(dtw_events)

        self.empty = False

    def calc_events(self):

        if self.read.has_events:
            self.dtw['sum'] = self.dtw['signal'] * self.dtw['length']

        grp = self.dtw.groupby("ref")
        sigs = grp['signal']

        ref_coords = np.abs(self.bcaln.y_min + grp['ref'].first())

        self.events = pd.DataFrame({
            "ref"   : ref_coords,
            "start"  : grp['sample'].min(),
            "kmer" : grp['kmer'].first(),
        })#.reset_index(drop=True)

        if self.read.has_events:
            self.events['length'] = grp['length'].sum()
            self.events['mean'] = grp['sum'].sum() / self.events['length']
            self.events['stdv'] = grp['signal'].std().fillna(0) #TODO not legit
            self.dtw.drop(columns=['sum'])
        else:
            self.events['length'] = grp['signal'].count()
            self.events['mean'] = grp['signal'].mean() 
            self.events['stdv'] = grp['signal'].std().fillna(0)

        #self.events.set_index('ref', inplace=True)
        #self.events.sort_index(inplace=True)

    def load_kmers(self):

        shift = nt.K - 1
        st = self.ref_min - (shift if not self.bcaln.flip_ref else 0)
        en = self.ref_max + (shift if self.bcaln.flip_ref else 0)

        #ref_len = self.bcaln.rf_en-self.bcaln.rf_st+K  
        #shift = K - 1                                  
        #if not self.bcaln.flip_ref:                    
        #    st = self.bcaln.rf_st - shift              
        #    en = st + ref_len                          
        #else:                                          
        #    en = self.bcaln.rf_en + shift              
        #    st = en - ref_len


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

        self.ref_kmers = np.insert(kmers, 0, [0]*pad)
        #self.ref_kmers = kmers
        #print([unc.kmer_to_str(k) for k in kmers])


    def get_dtw_args(self, read_block, ref_start, ref_kmers):
        common = (read_block['norm_sig'].to_numpy(), ref_kmers, self.model)
        qry_len = len(read_block)
        ref_len = len(ref_kmers)

        #TODO slow, probably bc of "searchsorted" line
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
                if r <= self.bcaln.df.loc[tgt,'ref'] - ref_start:
                    r += 1
                else:
                    q += 1

            return common + (self.prms.band_width, band_lls)

        elif self.method == "StaticBDTW":
            return common + (self.prms.band_width, self.prms.band_shift)

        else:
            return common + (DTW_GLOB,)

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


    def calc_dtw(self):
        self.mats = list()

        path_qrys = list()
        path_refs = list()

        band_blocks = list()

        block_min = self.bcaln.df['sample'].searchsorted(self.samp_min)
        block_max = self.bcaln.df['sample'].searchsorted(self.samp_max)


        y_min = self.bcaln.df['ref'][block_min]

        block_starts = np.insert(self.bcaln.ref_gaps, 0, block_min)
        block_ends   = np.append(self.bcaln.ref_gaps, block_max)

        for st, en in [(block_min, block_max)]:#zip(block_starts, block_ends):
        #for st, en in zip(block_starts, block_ends):
            samp_st = self.bcaln.df.loc[st,'sample']
            samp_en = self.bcaln.df.loc[en-1,'sample']

            ref_st = self.bcaln.df.loc[st,"ref"]
            ref_en = self.bcaln.df.loc[en-1,"ref"]

            read_block = self.read.sample_range(samp_st, samp_en)

            block_signal = read_block['norm_sig'].to_numpy()
            block_kmers = self.ref_kmers[ref_st-y_min:ref_en-y_min]

            args = self.get_dtw_args(read_block, ref_st, block_kmers)

            dtw = self.dtw_fn(*args)

            #TODO flip in traceback
            path = np.flip(dtw.path)
            path_qrys.append(read_block.index[path['qry']])
            path_refs.append(ref_st + path['ref'])

            if hasattr(dtw, "ll"):
                band_blocks.append(
                    self.ll_to_df(dtw.ll, read_block, ref_st, len(block_kmers))
                )

        self.dtw = pd.DataFrame({'ref': np.concatenate(path_refs)}, 
                               index = np.concatenate(path_qrys),
                               dtype='Int32') \
                  .join(self.read.df) \
                  .drop(columns=['mean', 'stdv', 'mask'], errors='ignore') \
                  .rename(columns={'start' : 'sample', 'norm_sig' : 'signal'})
        self.dtw['kmer'] = self.ref_kmers[self.dtw['ref'].astype(int)-y_min]

        #self.dtw['ref'] += self.bcaln.y_min

        if len(band_blocks) == 0:
            self.bands = None
        elif len(band_blocks) > 1:
            self.bands = pd.concat(band_blocks)
        else:
            self.bands = pd.DataFrame(band_blocks[0])

    def load_dtw_events(self, event_file):
        self.events = pd.read_pickle(event_file).reset_index()

        block_min = self.bcaln.df['sample'].searchsorted(self.samp_min)
        y_min1 = self.bcaln.df['ref'][block_min]

        y_min = self.events['ref'].min()
        y_max = self.events['ref'].max()

        #block_min2 = self.events['start'].searchsorted(self.samp_min)
        #print(self.events)
        #print(block_min2, self.samp_min)
        #y_min2 = self.events['ref'].iloc[block_min2]


        self.events = self.events.loc[(self.events['start'] >= self.samp_min) & (self.events['start'] <= self.samp_max)]#.reset_index(drop=True)

        y_min2 = self.events['ref'].min()

        if self.bcaln.flip_ref:
            self.events['idx'] = -self.events['ref'] + y_max
        else:
            self.events['idx'] = self.events['ref'] - y_min


        self.events.set_index('idx', inplace=True)
        self.events.sort_index(inplace=True)

        self.dtw = self.events.drop(columns=["ref"]) \
                              .reset_index() \
                              .rename(columns={'idx' : 'ref', 'start' : 'sample', 'mean' : 'signal'})
        self.dtw.reset_index()
        self.bands = None

    def plot_dotplot(self, ax):
        if self.bands is not None:
            ax.fill_between(self.bands['samp'], self.bands['ref_st']-1, self.bands['ref_en'], zorder=1, color='#ccffdd', linewidth=1, edgecolor='black', alpha=0.5)


        #return ax.scatter(self.dtw['sample'], self.dtw['ref'],s=7,color="purple", zorder=2)
        return ax.step(self.dtw['sample'], self.dtw['ref'],where="post",color="purple", zorder=3, linewidth=3)

    def plot_dtw_events(self, ax_sig, ax_padiff):
        c = 'purple'

        samps = np.arange(self.samp_min, self.samp_max)

        raw_norm = np.zeros(len(samps))
        i = self.read.norm_params["end"].searchsorted(self.samp_min)
        while i < len(self.read.norm_params):
            n = self.read.norm_params.iloc[i]

            st = int(n["start"])
            if st >= self.samp_max: break

            en = int(n["end"])-1

            st = max(st, self.samp_min)
            en = min(en, self.samp_max)

            raw = self.read.f5.signal[st:en]
            raw_norm[st-self.samp_min:en-self.samp_min] = (n["scale"] * raw) + n["shift"]
            i += 1

        ymin = np.min(raw_norm[raw_norm>0])
        ymax = np.max(raw_norm[raw_norm>0])
        bases = nt.kmer_base(self.dtw['kmer'], 2)

        samp_bases = np.zeros(len(samps), int)
        for i in range(len(self.dtw)):
            st = int(self.dtw.iloc[i]['sample']-self.samp_min)
            en = int(st + self.dtw.iloc[i]['length'])
            samp_bases[st:en] = bases[i]

            
        def plot_base(base, color):
            ax_sig.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color, interpolate=True)
        plot_base(0, "#80ff80")
        plot_base(1, "#8080ff")
        plot_base(2, "#ffbd00")
        plot_base(3, "#ff8080")

        ax_sig.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, alpha=0.75, c="#777777")

        model_means = self.model.get_mean(self.events['kmer'])

        ax_sig.step(self.events['start'], self.model.get_mean(self.events['kmer']), color='white', linewidth=2, where="post")

        ax_sig.vlines(self.events['start'], ymin, ymax, linewidth=2, color="white")

        evts = (self.read.df['start'] >= self.samp_min) & (self.read.df['start'] < self.samp_max) & (self.read.df['norm_sig'] > 0)

        if self.read.has_events:
            ax_sig.step(self.read.df['start'][evts], self.read.df['norm_sig'][evts], where='post', color='black', linewidth=3)
        else:
            ax_sig.scatter(self.read.df['start'][evts], self.read.df['norm_sig'][evts], s=5, alpha=0.75, c="#777777") #TODO really need to store constants, this is so badly organized

        model_means = self.model.get_mean(self.events['kmer'])

        #ax_padiff.scatter(self.dtw['signal'], self.dtw['ref'], color='purple', s=5, alpha=0.25)

        pa_diffs = np.abs(self.events['mean'] - self.model.get_mean(self.events['kmer']))

        #ax_padiff.plot(model_means, self.events.index, color='forestgreen')

        ax_padiff.step(pa_diffs, self.events.index, color=c, where="post")

        #replace .events['ref'] with index with .index to fix for cmd for some reason

#TODO move into GuidedDTW, or ReadAln?
def save(dtw, track):
    events_out = dtw.events.reset_index(drop=True).set_index('ref').sort_index()

    track.add_read(dtw.read.id, dtw.read.f5.filename, events_out)

class BcFast5Aln:
    BCE_K = 4
    CIG_OPS_STR = "MIDNSHP=X"
    CIG_RE = re.compile("(\d+)(["+CIG_OPS_STR+"])")
    CIG_OPS = set(CIG_OPS_STR)
    CIG_INCR_ALL = {'M','=', 'X'}
    CIG_INCR_RD = CIG_INCR_ALL | {'I','S'}
    CIG_INCR_RF = CIG_INCR_ALL | {'D','N'}

    SUB = 0
    INS = 1
    DEL = 2
    ERR_TYPES = [SUB, INS, DEL]
    ERR_MARKS = ['o', 'P', '_']
    ERR_SIZES = [100, 150, 150]
    ERR_WIDTHS = [0,0,5]

    def __init__(self, read, paf):
        self.seq_fwd = read.conf.read_buffer.seq_fwd

        self.refgap_bps = list()
        self.sub_bps = list()
        self.ins_bps = list()
        self.del_bps = list()
        self.err_bps = None

        self.empty = (
                paf is None or 
                not read.f5.bc_loaded or 
                (not self.parse_cs(paf) and
                 not self.parse_cigar(paf))
        )
        if self.empty: return

        self.rf_name = paf.rf_name
        self.rf_st = paf.rf_st
        self.rf_en = paf.rf_en
        self.is_fwd = paf.is_fwd

        #TODO make c++ build this 
        moves = np.array(read.f5.moves, bool)
        bce_qrs = np.cumsum(read.f5.moves)
        bce_samps = read.f5.template_start + np.arange(len(bce_qrs)) * read.f5.bce_stride

        samp_bps = pd.DataFrame({
            'sample' : bce_samps,#[moves],
            'bp'     : np.cumsum(read.f5.moves),#[moves],
        })

        self.df = samp_bps.join(self.bp_rf_aln, on='bp').dropna()
        self.df.reset_index(inplace=True, drop=True)

        #Make sure ref starts at 0
        self.df['ref'] -= self.df['ref'].min()
        
        if self.err_bps is not None:
            self.errs = samp_bps.join(self.err_bps.set_index('bp'), on='bp').dropna()
            self.errs.reset_index(inplace=True, drop=True)
        else:
            self.errs = None

        self.ref_gaps = self.df[self.df['bp'].isin(self.refgap_bps)].index

        self.subs = self.df[self.df['bp'].isin(self.sub_bps)].index
        self.inss = self.df[self.df['bp'].isin(self.ins_bps)].index
        self.dels = self.df[self.df['bp'].isin(self.del_bps)].index

        #print((self.bce_samps == self.df[self.df['ref'] == self.bce_refs]['sample']).all())
        #self.rfgap_bces = self.bce_refs[self.bce_bps == 

        self.empty = len(self.df) == 0
        if self.empty: return

        self.flip_ref = paf.is_fwd != self.seq_fwd
        if self.flip_ref:
            rf_st = paf.rf_en - self.df['ref'].max() - 1
            rf_en = paf.rf_en
        else:
            rf_st = paf.rf_st 
            rf_en = paf.rf_st + self.df['ref'].max() + 1

        self.ref_bounds = (
            paf.rf_name,
            rf_st,
            rf_en,
            paf.is_fwd
        )

        self.y_min = -paf.rf_en if self.flip_ref else paf.rf_st
        self.y_max = self.y_min + self.df['ref'].max()

    def ref_tick_fmt(self, ref, pos=None):
        return int(np.round(np.abs(self.y_min + ref)))

    def parse_cs(self, paf):
        cs = paf.tags.get('cs', (None,)*2)[0]
        if cs is None: return False

        sys.stderr.write("Loading cs tag\n")

        #TODO rename to general cig/cs
        bp_rf_aln = list()
        err_bps = list()

        if self.seq_fwd:
            qr_i = paf.qr_st
            #rf_i = paf.rf_st
        else:
            qr_i = paf.qr_len - paf.qr_en 
            #rf_i = -paf.rf_en+1
        rf_i = 0

        cs_ops = re.findall("(=|:|\*|\+|-|~)([A-Za-z0-9]+)", cs)

        if paf.is_fwd != self.seq_fwd:
            cs_ops = reversed(cs_ops)

        for op in cs_ops:
            c = op[0]
            if c in {'=',':'}:
                l = len(op[1]) if c == '=' else int(op[1])
                bp_rf_aln += zip(range(qr_i, qr_i+l), range(rf_i, rf_i+l))
                qr_i += l
                rf_i += l

            elif c == '*':
                self.sub_bps.append(qr_i)
                bp_rf_aln.append((qr_i,rf_i))
                err_bps.append( (qr_i,rf_i,self.SUB) )
                qr_i += 1
                rf_i += 1

            elif c == '-':
                self.ins_bps.append(qr_i)
                err_bps.append( (qr_i,rf_i,self.DEL) )
                l = len(op[1])
                rf_i += l

            elif c == '+':
                self.del_bps.append(qr_i)
                err_bps.append( (qr_i,rf_i,self.INS) )

                l = len(op[1])
                qr_i += l

            elif c == '~':
                l = int(op[1][2:-2])
                self.refgap_bps.append(qr_i)
                rf_i += l

            else:
                print("UNIMPLEMENTED ", op)

        self.bp_rf_aln = pd.DataFrame(bp_rf_aln, columns=["bp","ref"], dtype='Int64')
        self.bp_rf_aln.set_index("bp", inplace=True)
        #TODO type shouldn't have to be 64 bit
        self.err_bps = pd.DataFrame(err_bps, columns=["bp","ref","type"], dtype='Int64')

        return True        


    def parse_cigar(self, paf):
        cig = paf.tags.get('cg', (None,)*2)[0]
        if cig is None: return False

        bp_rf_aln = list()#defaultdict(list)
        self.refgap_bps = list()

        if self.seq_fwd:
            qr_i = paf.qr_st
            #rf_i = paf.rf_st
        else:
            qr_i = paf.qr_len - paf.qr_en 
            #rf_i = -paf.rf_en+1

        rf_i = 0

        cig_ops = self.CIG_RE.findall(cig)

        if paf.is_fwd != self.seq_fwd:
            cig_ops = list(reversed(cig_ops))

        for l,c in cig_ops:
            l = int(l)
            incr_qr = c in self.CIG_INCR_RD
            incr_rf = c in self.CIG_INCR_RF
            qr_j = qr_i + (l if incr_qr else 1)
            rf_j = rf_i + (l if incr_rf else 1)

            if c == "M":
                bp_rf_aln += zip(range(qr_i, qr_j), range(rf_i, rf_j))
            elif c == "N":
                self.refgap_bps.append(qr_i)

            if incr_qr:
                qr_i = qr_j 

            if incr_rf:
                rf_i = rf_j 

        self.bp_rf_aln = pd.DataFrame(bp_rf_aln, columns=["bp","ref"], dtype='Int64')
        self.bp_rf_aln.set_index("bp", inplace=True)

        return True

    def get_xy(self, i):
        df = self.df.loc[i]
        return (df['sample'], df['ref']-0.5)

    def plot_scatter(self, ax, real_start=False, samp_min=None, samp_max=None):
        if samp_min is None: samp_min = 0
        if samp_max is None: samp_max = self.df['sample'].max()
        i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        return ax.scatter(self.df['sample'][i], self.df['ref'][i], color='orange', zorder=2,s=20)

    def plot_step(self, ax, real_start=False, samp_min=None, samp_max=None):
        i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        ret = ax.step(self.df['sample'][i], self.df['ref'][i], color='orange', zorder=1, where='post')

        if self.errs is not None:
            for t in self.ERR_TYPES:
                e = self.errs[self.errs['type'] == t]
                ax.scatter(
                    e['sample'], e['ref'], 
                    color='red', zorder=3, 
                    s=self.ERR_SIZES[t],
                    linewidth=self.ERR_WIDTHS[t],
                    marker=self.ERR_MARKS[t]
                )

        return ret

    def samp_to_ref(self, samp):
        if self.flip_ref: ref = -ref
        ref = ref - self.y_min
        i = np.clip(self.df['ref'].searchsorted(ref), 0, len(self.df)-1)
        return self.df['sample'][i]

    def ref_to_samp(self, ref):
        if self.flip_ref: ref = -ref
        ref = ref - self.y_min
        i = np.clip(self.df['ref'].searchsorted(ref), 0, len(self.df)-1)
        return self.df['sample'][i]

class Fast5Processor(Fast5Reader):
    def __next__(self):
        return ProcRead(Fast5Reader.__next__(self), conf=self.conf)
    
    def __getitem__(self, read_id):
        return ProcRead(Fast5Reader.__getitem__(self, read_id), conf=self.conf)
