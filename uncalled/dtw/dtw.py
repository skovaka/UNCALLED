#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict
import re
import time
from typing import NamedTuple
from matplotlib.colors import Normalize
import pandas as pd

from ..pafstats import parse_paf, PafEntry
from ..config import Config
from _uncalled import PORE_MODELS

class ReadAln:

    def __init__(self, index, aln, is_rna=False, ref_bounds=None):
        if type(aln) != PafEntry:
            raise RuntimeError("ReadAlns can only be initialized from PafEntrys currently")

        self.index = index
        self.is_rna = is_rna

        self.clipped = False

        if ref_bounds is None:
            self.ref_bounds = (aln.rf_name, aln.rf_st, aln.rf_en, aln.is_fwd)
        else:
            if aln.rf_st < ref_bounds[1]:
                ref_st = ref_bounds[1]
                clipped = True
            else:
                ref_st = aln.rf_st

            if aln.rf_en > ref_bounds[2]:
                ref_en = ref_bounds[2]
                clipped = True
            else:
                ref_en = aln.rf_en

            if ref_st > ref_en:
                self.empty = True
                return

            self.ref_bounds = (aln.rf_name, ref_st, ref_en, aln.is_fwd)

        self._init_mirror_coords()

        self.empty = False

    def miref_to_ref(self, miref):
        return index.mirror_to_ref(miref) 

    def ref_to_miref(self, ref):
        return self.index.mirror_ref_coords(self.ref_name, ref, ref, self.is_fwd, self.is_rna)[0]

    def _init_mirror_coords(self):
        self.miref_start, self.miref_end = self.index.mirror_ref_coords(*self.ref_bounds, self.is_rna)


    @property
    def ref_name(self):
        """The sequence name of the alignment reference coordinate"""
        return self.ref_bounds[0]

    @property
    def ref_start(self):
        """The start of the alignment reference coordinate"""
        return self.ref_bounds[1]

    @property
    def ref_end(self):
        """The end of the alignment reference coordinate"""
        return self.ref_bounds[2]

    @property
    def is_fwd(self):
        return self.ref_bounds[3]

    def get_index_kmers(self, index, kmer_shift=4):
        """Returns the k-mer sequence at the alignment reference coordinates"""
        start = self.miref_start - kmer_shift

        if start < 0:
            lpad = -self.miref_start
            start = 0
        else:
            lpad = 0

        kmers = index.get_kmers(start, self.miref_end, self.is_rna)

        return np.insert(kmers, 0, [0]*lpad)

    
class Track:
    CONF_FNAME = "conf.toml"
    ALN_DIR = "alns"
    ALN_SUFFIX = ".pkl"

    INDEX_FNAME = "filename_mapping.txt"
    INDEX_HEADER = "read_id\tfilename\taln_file"

    #REF_LOCS = "read_id\tsample_start\tsample_end\tref_name\tref_start\tref_end\tstrand"

    #TODO make fast5_file alias for filename in Fast5Reader
    #also, in fast5 reader make static fast5 filename parser

    WRITE_MODE = "w"
    READ_MODE = "r"
    MODES = {WRITE_MODE, READ_MODE}

    def __init__(self, path, mode="r", conf=None, overwrite=False):
        self.path = path.strip("/")
        self.mode = mode

        self.conf = conf if conf is not None else Config()

        if mode == self.WRITE_MODE:
            os.makedirs(self.aln_dir, exist_ok=overwrite)

        self.index_file = open(self.index_filename, mode)

        if mode == self.READ_MODE:
            self.conf.load_toml(self.config_filename)
            self.conf.fast5_reader.fast5_index = self.index_filename
            self._load_index()
            #self._load_reads()

        elif mode == self.WRITE_MODE:
            self.conf.to_toml(self.config_filename)
            self.index_file.write(self.INDEX_HEADER + "\n")

    @property
    def read_ids(self):
        return list(self.index.index)

    def add_read(self, read_id, fast5_fname, rae_df):
        if self.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")

        aln_fname = self.aln_fname(read_id)
        rae_df.to_pickle(aln_fname)

        self.index_file.write("\t".join([read_id, fast5_fname, aln_fname]) + "\n")

    def _load_index(self):
        self.index = pd.read_csv(self.index_file, sep="\t", index_col="read_id")

    def get_aln(self, read_id, ref_bounds=None):
        aln = pd.read_pickle(self.aln_fname(read_id))

        if ref_bounds is not None:
            _,st,en = ref_bounds[:3]
            return aln[st:en]
        return aln

    def get_matrix(self, ref_bounds, mm2_paf=None, partial_overlap=False):

        if mm2_paf is None:
            mm2_paf = self.conf.align.mm2_paf

        read_filter = set(self.index.index)
        mm2s = {p.qr_name : p
                 for p in parse_paf(
                    mm2_paf,
                    ref_bounds,
                    read_filter=read_filter,
                    full_overlap=not partial_overlap
        )}

        mat = TrackMatrix(self, ref_bounds, mm2s, conf=self.conf)

        for read_id,read in self.index.iterrows():
            mm2 = mm2s.get(read_id, None)
            if mm2 is None: continue
            df = pd.read_pickle(read["aln_file"])
            mat._add_read(df, mm2)

        mat._flatten()

        return mat

    def close(self):
        self.index_file.close()

    @property
    def config_filename(self):
        return os.path.join(self.path, self.CONF_FNAME)

    @property
    def index_filename(self):
        return os.path.join(self.path, self.INDEX_FNAME)
    
    @property
    def aln_dir(self):
        return os.path.join(self.path, self.ALN_DIR)
    
    def aln_fname(self, read_id):
        return os.path.join(self.aln_dir, read_id+self.ALN_SUFFIX)

class TrackMatrix:

    KMER_LAYER = 0
    PA_LAYER = 1
    DWELL_LAYER = 2
    PA_DIFF_LAYER = 3

    def __init__(self, track, ref_bounds, mm2s, height=None, conf=None):
        self.conf = conf if conf is not None else Config()
        self.track = track

        self.ref_bounds = ref_bounds
        self.width = self.ref_end-self.ref_start
        self.height = height

        self._layers = defaultdict(list)
        self.reads = defaultdict(list)
        self.mask = list()
        self.mm2s = dict()

        model_name = self.conf.mapper.pore_model
        #TODO probably need to rethink fwd/rev compl, but either way clean this up
        if model_name.endswith("_compl"):
            model_name = model_name[:-5]+"templ"

        self.model = PORE_MODELS[model_name]

        self.ref_to_x = pd.Series(
                np.arange(self.width, dtype=int),
                index=pd.RangeIndex(self.ref_start, self.ref_end)
        )

    def __getitem__(self, i):
        return self._layers.__getitem__(i)
    
    @property
    def kmer(self):
        return self._layers[self.KMER_LAYER]

    @property
    def pa(self):
        return self._layers[self.PA_LAYER]
    
    @property
    def dwell(self):
        return self._layers[self.DWELL_LAYER]
    
    @property
    def pa_diff(self):
        return self._layers[self.PA_DIFF_LAYER]

    def _add_read(self, df, mm2_paf):

        dtw_roi = df.loc[self.ref_start:self.ref_end-1]

        xs = self.ref_to_x.reindex(self.ref_to_x.index.intersection(dtw_roi.index))

        roi_mask = np.zeros(self.width)
        roi_mask[xs] = True
        self.mask.append(roi_mask)

        pa_diffs = self.model.match_diff(dtw_roi['mean'], dtw_roi['kmer'])
        dwell = 1000 * dtw_roi['length'] / self.conf.read_buffer.sample_rate

        self._add_layer_row(self.KMER_LAYER, dtw_roi['kmer'], xs)
        self._add_layer_row(self.PA_LAYER, dtw_roi['mean'], xs)
        self._add_layer_row(self.PA_DIFF_LAYER, pa_diffs, xs)
        self._add_layer_row(self.DWELL_LAYER, dwell, xs)

        self.reads['ref_start'].append(df.index.min())
        self.reads['id'].append(mm2_paf.qr_name)
        self.reads['fwd'].append(mm2_paf.is_fwd)

        self.mm2s[mm2_paf.qr_name] = mm2_paf

    def _add_layer_row(self, layer, vals, xs):
        row = np.zeros(self.width)
        row[xs] = vals
        self._layers[layer].append(row)
        return row

    def _flatten(self):
        self.reads = pd.DataFrame(self.reads) \
                     .sort_values(['fwd', 'ref_start'], ascending=[False, True])

        self.has_fwd = np.any(self.reads['fwd'])
        self.has_rev = not np.all(self.reads['fwd'])

        read_order = self.reads.index.to_numpy()

        layer_names = sorted(self._layers.keys())

        self._layers = np.stack([
            np.stack(self._layers[l])[read_order] for l in layer_names
        ])
        self.mask = np.stack(self.mask)[read_order].astype(bool)

        self.height = len(self.reads)

        self.norms = [Normalize(np.min(l[self.mask]), np.max(l[self.mask])) for l in self._layers]
        for l in [self.PA_LAYER, self.DWELL_LAYER]:
            lmask = self._layers[l][self.mask]
            self.norms[l].vmax = min(
                lmask.max(),
                np.median(lmask) + 2 * lmask.std()
            )

    def sort(self, layer, ref):
        order = np.argsort(self._layers[layer,:,ref])
        self._layers = self._layers[:,order,:]
        self.reads = self.reads.iloc[order]

    @property
    def ref_name(self):
        return self.ref_bounds[0]

    @property
    def ref_start(self):
        return self.ref_bounds[1]

    @property
    def ref_end(self):
        return self.ref_bounds[2]

def ref_coords(coord_str):
    print("COORD", coord_str)
    spl = coord_str.split(":")
    ch = spl[0]
    st,en = spl[1].split("-")

    coord = (ch, int(st), int(en))

    if len(spl) == 2:
        return coord
    else:
        return coord + (spl[2] == "+",)

class BcFast5Aln(ReadAln):
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

    def __init__(self, index, read, paf, ref_bounds=None):
        self.seq_fwd = read.conf.read_buffer.seq_fwd
        ReadAln.__init__(self, index, paf, not self.seq_fwd, ref_bounds=ref_bounds)

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

        #TODO make c++ build this 
        moves = np.array(read.f5.moves, bool)
        bce_qrs = np.cumsum(read.f5.moves)
        bce_samps = read.f5.template_start + np.arange(len(bce_qrs)) * read.f5.bce_stride

        samp_bps = pd.DataFrame({
            'sample' : bce_samps,#[moves],
            'bp'     : np.cumsum(read.f5.moves),#[moves],
        })

        self.df = samp_bps.join(self.bp_miref_aln, on='bp').dropna()
        self.df.reset_index(inplace=True, drop=True)

        if self.err_bps is not None:
            self.errs = samp_bps.join(self.err_bps.set_index('bp'), on='bp').dropna()
            self.errs.reset_index(inplace=True, drop=True)
        else:
            self.errs = None

        self.ref_gaps = self.df[self.df['bp'].isin(self.refgap_bps)].index

        self.subs = self.df[self.df['bp'].isin(self.sub_bps)].index
        self.inss = self.df[self.df['bp'].isin(self.ins_bps)].index
        self.dels = self.df[self.df['bp'].isin(self.del_bps)].index

        self.empty = len(self.df) == 0
        if self.empty: return

        self.flip_ref = paf.is_fwd != self.seq_fwd

        #if self.flip_ref:
        #    ref_st = paf.rf_en - self.df['miref'].max() - 1
        #    ref_en = paf.rf_en
        #else:
        #    ref_st = paf.rf_st 
        #    ref_en = paf.rf_st + self.df['miref'].max() + 1

        #self.ref_bounds = (
        #    paf.rf_name,
        #    ref_st,
        #    ref_en,
        #    paf.is_fwd
        #)

        self.y_min = -paf.rf_en if self.flip_ref else paf.rf_st
        self.y_max = self.y_min + self.df['miref'].max()

    def ref_tick_fmt(self, ref, pos=None):
        return int(np.round(np.abs(self.y_min + ref)))

    def parse_cs(self, paf):
        cs = paf.tags.get('cs', (None,)*2)[0]
        if cs is None: return False

        sys.stderr.write("Loading cs tag\n")

        #TODO rename to general cig/cs
        bp_miref_aln = list()
        err_bps = list()

        if self.seq_fwd:
            qr_i = paf.qr_st
            #rf_i = paf.rf_st
        else:
            qr_i = paf.qr_len - paf.qr_en 
            #rf_i = -paf.rf_en+1

        print(self.miref_start, self.clipped, "MIORR")
        mr_i = self.miref_start

        cs_ops = re.findall("(=|:|\*|\+|-|~)([A-Za-z0-9]+)", cs)

        if paf.is_fwd != self.seq_fwd:
            cs_ops = reversed(cs_ops)

        for op in cs_ops:
            c = op[0]
            if c in {'=',':'}:
                l = len(op[1]) if c == '=' else int(op[1])
                bp_miref_aln += zip(range(qr_i, qr_i+l), range(mr_i, mr_i+l))
                qr_i += l
                mr_i += l

            elif c == '*':
                self.sub_bps.append(qr_i)
                bp_miref_aln.append((qr_i,mr_i))
                err_bps.append( (qr_i,mr_i,self.SUB) )
                qr_i += 1
                mr_i += 1

            elif c == '-':
                self.ins_bps.append(qr_i)
                err_bps.append( (qr_i,mr_i,self.DEL) )
                l = len(op[1])
                mr_i += l

            elif c == '+':
                self.del_bps.append(qr_i)
                err_bps.append( (qr_i,mr_i,self.INS) )

                l = len(op[1])
                qr_i += l

            elif c == '~':
                l = int(op[1][2:-2])
                self.refgap_bps.append(qr_i)
                mr_i += l

            else:
                print("UNIMPLEMENTED ", op)

        self.bp_miref_aln = pd.DataFrame(bp_miref_aln, columns=["bp","miref"], dtype='Int64')
        self.bp_miref_aln.set_index("bp", inplace=True)

        #TODO type shouldn't have to be 64 bit
        self.err_bps = pd.DataFrame(err_bps, columns=["bp","miref","type"], dtype='Int64')

        return True        


    def parse_cigar(self, paf):
        cig = paf.tags.get('cg', (None,)*2)[0]
        if cig is None: return False

        bp_miref_aln = list()#defaultdict(list)
        self.refgap_bps = list()

        if self.seq_fwd:
            qr_i = paf.qr_st
            mr_i = self.ref_to_miref(paf.rf_st)
        else:
            qr_i = paf.qr_len - paf.qr_en 
            mr_i = self.ref_to_miref(paf.rf_en)

        #mr_i = self.miref_start
        print(self.miref_start, self.clipped, "MIORR")
        mr_bounds = range(self.miref_start, self.miref_end)

        cig_ops = self.CIG_RE.findall(cig)

        if paf.is_fwd != self.seq_fwd:
            cig_ops = list(reversed(cig_ops))

        for l,c in cig_ops:
            l = int(l)
            incr_qr = c in self.CIG_INCR_RD
            incr_rf = c in self.CIG_INCR_RF
            qr_j = qr_i + (l if incr_qr else 1)
            mr_j = mr_i + (l if incr_rf else 1)

            if c == "M":
                for qr, mr in zip(range(qr_i, qr_j), range(mr_i, mr_j)):
                    if mr in mr_bounds:
                        bp_miref_aln.append((qr,mr))
                #bp_miref_aln += zip(range(qr_i, qr_j), range(mr_i, mr_j))
            elif c == "N":
                if mr_i in mr_bounds:
                    bp_miref_aln.append((qr_i,mr))

            if incr_qr:
                qr_i = qr_j 

            if incr_rf:
                mr_i = mr_j 

        self.bp_miref_aln = pd.DataFrame(bp_miref_aln, columns=["bp","miref"], dtype='Int64')
        self.bp_miref_aln.set_index("bp", inplace=True)

        return True

    def get_xy(self, i):
        df = self.df.loc[i]
        return (df['sample'], df['miref']-0.5)

    def plot_scatter(self, ax, real_start=False, samp_min=None, samp_max=None):
        if samp_min is None: samp_min = 0
        if samp_max is None: samp_max = self.df['sample'].max()
        i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        return ax.scatter(self.df['sample'][i], self.df['miref'][i], color='orange', zorder=2,s=20)

    def plot_step(self, ax, real_start=False, samp_min=None, samp_max=None):
        i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        ret = ax.step(self.df['sample'][i], self.df['miref'][i], color='orange', zorder=1, where='post')

        if self.errs is not None:
            for t in self.ERR_TYPES:
                e = self.errs[self.errs['type'] == t]
                ax.scatter(
                    e['sample'], e['miref'], 
                    color='red', zorder=3, 
                    s=self.ERR_SIZES[t],
                    linewidth=self.ERR_WIDTHS[t],
                    marker=self.ERR_MARKS[t]
                )

        return ret

    def ref_to_samp(self, ref):
        return self.miref_to_samp(self.ref_to_miref(ref))
        

    def miref_to_samp(self, miref):
        i = np.clip(self.df['miref'].searchsorted(miref), 0, len(self.df)-1)
        return self.df['sample'][i]
