#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict, namedtuple
import re
import time
from matplotlib.colors import Normalize
import pandas as pd
import scipy.stats
import copy

import matplotlib.pyplot as plt

from ..pafstats import parse_paf, PafEntry
from ..config import Config, Opt
from .. import BwaIndex, nt, PoreModel

class RefCoord:
    def __init__(self, name, start=None, end=None, fwd=None):
        self.fwd = fwd
        if start is None and end is None:
            if isinstance(name ,str):
                self._init_str(name)
            elif isinstance(name, tuple):
                self._init_tuple(name)

        elif start is None:
            raise ValueError("RefCoords must include a start coordinate")
        
        else:
            self.name = name
            self.start = start
            self.end = end

    def _init_str(self, coord_str):
        spl = coord_str.split(":")
        self.name = spl[0]
        coords = tuple(map(int, spl[1].split("-")))

        if len(coords) == 2:
            self.start, self.end = coords
        elif len(coords) == 1:
            self.start, = coords
            self.end = None
        else:
            raise ValueError("RefCoords must contain one or two coordinate")

        if len(spl) == 3:
            self.fwd = spl[2] == "+"
        else:
            self.fwd = None

    def _init_tuple(self, coords):
        self.name = coords[0]
        self.start = coords[1]
        if len(coords) > 2:
            self.end = coords[2]
        if len(coords) > 3:
            self.end = coords[3]

    def __repr__(self):
        s = "%s:%d" % (self.name, self.start)
        if self.end is not None:
            s += "-%d" % self.end
        if self.fwd is not None:
            s += " (%s)" % ("+" if self.fwd else "-")
        return s

LayerMeta = namedtuple("LayerMeta", ["type", "label"])

LAYER_META = {
    "ref"     : LayerMeta(int, "Reference Coordinate"),
    "start"   : LayerMeta(int, "Sample Start"),
    "length"  : LayerMeta(int, "Sample Length"),
    "current" : LayerMeta(float, "Mean Current (pA)"),
    "kmer"    : LayerMeta(int, "Reference K-mer"),
    "refmir"  : LayerMeta(int, "Mirrored Packed Ref. Coord."),
}

class ReadAln:

    def __init__(self, index, aln, df=None, is_rna=False, ref_bounds=None):
        if not isinstance(aln, PafEntry):
            raise RuntimeError("ReadAlns can only be initialized from PafEntrys currently")

        self.index = index
        self.read_id = aln.qr_name
        self.is_rna = is_rna

        self.clipped = False
    
        self.set_ref_bounds(aln, ref_bounds)

        if self.empty: 
            return

        self._init_mirror_coords()

        self.dfs = set()

        if df is not None:
            self.set_aln(df)

            has_ref = self.aln.index.name == "ref"
            has_refmir = "refmir" in self.aln.columns

            #TODO check for required columns
            if not has_ref and not has_refmir:
                raise RuntimeError("ReadAln DataFrame must include a column named \"%s\" or \"%s\"" % ("ref", "refmir"))
            
            if has_ref and not has_refmir:
                self.aln["refmir"] = self.index.ref_to_refmir(self.ref_id, self.aln.index, self.is_fwd, self.is_rna)

            elif not has_ref and has_refmir:
                self.aln[REF_COL] = self.index.mirref_to_ref(self.aln["refmir"])

            self.aln.sort_values("refmir", inplace=True)
        

    def set_ref_bounds(self, aln, ref_bounds):
        if ref_bounds is None:
            self.ref_bounds = RefCoord(aln.rf_name, aln.rf_st, aln.rf_en, aln.is_fwd)
        else:
            if aln.rf_st < ref_bounds.start:
                ref_st = ref_bounds.start
                clipped = True
            else:
                ref_st = aln.rf_st

            if aln.rf_en > ref_bounds.end:
                ref_en = ref_bounds.end
                clipped = True
            else:
                ref_en = aln.rf_en

            if ref_st > ref_en:
                self.empty = True
                return

            self.ref_bounds = RefCoord(aln.rf_name, ref_st, ref_en, aln.is_fwd)

        self.ref_id = self.index.get_ref_id(self.ref_name)

        self.empty = False

    #TODO ReadAln stores RefCoords, handles all this conversion?

    def refmir_to_ref(self, refmir):
        return self.index.refmir_to_ref(refmir) 

    def ref_to_refmir(self, ref):
        return self.index.ref_to_refmir(self.ref_name, ref, ref, self.is_fwd, self.is_rna)[0]

    def ref_to_samp(self, ref):
        return self.refmir_to_samp(self.ref_to_refmir(ref))
        
    def refmir_to_samp(self, refmir):
        i = np.clip(self.aln['refmir'].searchsorted(refmir), 0, len(self.aln)-1)
        return self.aln['sample'][i]
    
    def calc_refmir(self):
        self.aln["refmir"] = self.index.ref_to_refmir(self.ref_id, self.aln.index, self.is_fwd, self.is_rna)

    def _init_mirror_coords(self):
        #self.refmir_start, self.refmir_end = self.index.ref_to_refmir(*self.ref_bounds, self.is_rna)
        self.refmir_start, self.refmir_end = self.index.ref_to_refmir(self.ref_name, self.ref_start, self.ref_end, self.is_fwd, self.is_rna)

    @property
    def ref_name(self):
        """The sequence name of the alignment reference coordinate"""
        return self.ref_bounds.name

    @property
    def ref_start(self):
        """The start of the alignment reference coordinate"""
        return self.ref_bounds.start

    @property
    def ref_end(self):
        """The end of the alignment reference coordinate"""
        return self.ref_bounds.end

    @property
    def is_fwd(self):
        return self.ref_bounds.fwd
    
    def sort_ref(self):
        self.aln.sort_values("ref", inplace=True)

    def sort_refmir(self):
        self.aln.sort_values("refmir", inplace=True)

    def get_samp_bounds(self):
        samp_min = self.aln['start'].min()
        max_i = self.aln['start'].argmax()
        samp_max = self.aln['start'].iloc[max_i] + self.aln['length'].iloc[max_i]
        return samp_min, samp_max
    
    #def set_bands(self, bands):
    def set_df(self, df, name):
        self.dfs.add(name)
        df = df.loc[self.ref_start:self.ref_end-1].copy()
        setattr(self, name, df)


    def set_aln(self, df):
        self.set_df(df, "aln")

    def set_bcerr(self, df):
        self.set_df(df, "bcerr")

    def set_subevent_aln(self, aln, ref_mirrored=False, kmer_str=False, ref_col="refmir", start_col="start", length_col="length", mean_col="current", kmer_col="kmer"):

        aln["cuml_mean"] = aln[length_col] * aln[mean_col]

        grp = aln.groupby(ref_col)

        if kmer_str:
            kmers = [nt.kmer_rev(nt.str_to_kmer(k,0)) for k in grp[kmer_col].first()]
        else:
            kmers = grp[kmer_col].first()

        if ref_mirrored:
            refmirs = grp[ref_col].first()
            refs = self.refmir_to_ref(refmirs)
        else:
            refs = grp[ref_col].first()

        lengths = grp[length_col].sum()

        aln = pd.DataFrame({
            "ref"    : refs.astype("int64"),
            "kmer"   : kmers.astype("uint16"),
            "start"  : grp[start_col].min().astype("uint32"),
            "length" : lengths.astype("uint32"),
            "current"   : grp["cuml_mean"].sum() / lengths
        })
        
        if ref_mirrored:
            aln["refmir"] = refmirs

        self.set_aln(aln.set_index("ref").sort_values("ref"))

    def get_index_kmers(self, index, kmer_shift=4):
        """Returns the k-mer sequence at the alignment reference coordinates"""
        start = self.refmir_start - kmer_shift

        if start < 0:
            lpad = -start
            start = 0
        else:
            lpad = 0

        kmers = index.get_kmers(start, self.refmir_end, self.is_rna)
        return np.array(kmers)

        #return np.insert(kmers, 0, [0]*lpad)

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

    #Error = pd.Caegorical(["SUB", "INS", "DEL"])

    def __init__(self, index, read, paf, ref_bounds=None):
        self.seq_fwd = read.conf.read_buffer.seq_fwd #TODO just store is_rna
        ReadAln.__init__(self, index, paf, is_rna=not self.seq_fwd, ref_bounds=ref_bounds)
        if self.empty: return

        self.refgap_bps = list()
        self.sub_bps = list()
        self.ins_bps = list()
        self.del_bps = list()
        self.err_bps = None

        self.flip_ref = paf.is_fwd != self.seq_fwd

        self.empty = (
                not read.f5.bc_loaded or 
                (not self.parse_cs(paf) and
                 not self.parse_cigar(paf))
        )
        if self.empty: 
            return

        #TODO make c++ build this 
        moves = np.array(read.f5.moves, bool)
        bce_qrs = np.cumsum(read.f5.moves)
        bce_samps = read.f5.template_start + np.arange(len(bce_qrs)) * read.f5.bce_stride

        samp_bps = pd.DataFrame({
            'sample' : bce_samps,#[moves],
            'bp'     : np.cumsum(read.f5.moves),#[moves],
        })

        df = samp_bps.join(self.bp_refmir_aln, on='bp').dropna()
        df["ref"] = self.refmir_to_ref(df["refmir"])
        df.set_index("ref", drop=True, inplace=True)
        self.set_aln(df.sort_index())
        #self.aln.reset_index(inplace=True, drop=True)

        if self.err_bps is not None:
            self.errs = samp_bps.join(self.err_bps.set_index('bp'), on='bp').dropna()
            self.errs.reset_index(inplace=True, drop=True)
            self.errs["ref"] = self.refmir_to_ref(self.errs["refmir"])
            #self.errs.set_index("ref", inplace=True)
        else:
            self.errs = None


        self.ref_gaps = self.aln[self.aln['bp'].isin(self.refgap_bps)].index

        self.subs = self.aln[self.aln['bp'].isin(self.sub_bps)].index
        self.inss = self.aln[self.aln['bp'].isin(self.ins_bps)].index
        self.dels = self.aln[self.aln['bp'].isin(self.del_bps)].index

        self.empty = len(self.aln) == 0
        if self.empty: 
            return

        self.y_min = -paf.rf_en if self.flip_ref else paf.rf_st
        self.y_max = self.y_min + self.aln['refmir'].max()

    def parse_cs(self, paf):
        cs = paf.tags.get('cs', (None,)*2)[0]
        if cs is None: return False

        #TODO rename to general cig/cs
        bp_refmir_aln = list()
        err_bps = list()

        if self.seq_fwd:
            qr_i = paf.qr_st
            #rf_i = paf.rf_st
        else:
            qr_i = paf.qr_len - paf.qr_en 
            #rf_i = -paf.rf_en+1

        mr_i = self.refmir_start

        cs_ops = re.findall("(=|:|\*|\+|-|~)([A-Za-z0-9]+)", cs)

        if paf.is_fwd != self.seq_fwd:
            cs_ops = reversed(cs_ops)

        for op in cs_ops:
            c = op[0]
            if c in {'=',':'}:
                l = len(op[1]) if c == '=' else int(op[1])
                bp_refmir_aln += zip(range(qr_i, qr_i+l), range(mr_i, mr_i+l))
                qr_i += l
                mr_i += l

            elif c == '*':
                self.sub_bps.append(qr_i)
                bp_refmir_aln.append((qr_i,mr_i))
                err_bps.append( (qr_i,mr_i,"SUB",op[1][1].upper()) )
                qr_i += 1
                mr_i += 1

            elif c == '-':
                self.ins_bps.append(qr_i)
                err_bps.append( (qr_i,mr_i,"DEL",op[1].upper()) )
                l = len(op[1])
                mr_i += l

            elif c == '+':
                self.del_bps.append(qr_i)
                err_bps.append( (qr_i,mr_i,"INS",op[1].upper()) )

                l = len(op[1])
                qr_i += l

            elif c == '~':
                l = int(op[1][2:-2])
                self.refgap_bps.append(qr_i)
                mr_i += l

            else:
                print("UNIMPLEMENTED ", op)

        self.bp_refmir_aln = pd.DataFrame(bp_refmir_aln, columns=["bp","refmir"], dtype='Int64')
        self.bp_refmir_aln.set_index("bp", inplace=True)

        #TODO type shouldn't have to be 64 bit
        self.err_bps = pd.DataFrame(err_bps, columns=["bp","refmir","type","seq"])#, dtype='Int64')

        return True        

    def parse_cigar(self, paf):
        cig = paf.tags.get('cg', (None,)*2)[0]
        if cig is None: return False

        bp_refmir_aln = list()#defaultdict(list)
        self.refgap_bps = list()

        #mr_i = self.refmir_start
        if self.seq_fwd:
            qr_i = paf.qr_st
        else:
            qr_i = paf.qr_len - paf.qr_en 

        if self.flip_ref:
            mr_i = self.ref_to_refmir(paf.rf_en)
        else:
            mr_i = self.ref_to_refmir(paf.rf_st)

        mr_bounds = range(self.refmir_start, self.refmir_end)

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
                        bp_refmir_aln.append((qr,mr))
                #bp_refmir_aln += zip(range(qr_i, qr_j), range(mr_i, mr_j))
            elif c == "N":
                if mr_i in mr_bounds:
                    bp_refmir_aln.append((qr_i,mr))

            if incr_qr:
                qr_i = qr_j 

            if incr_rf:
                mr_i = mr_j 

        self.bp_refmir_aln = pd.DataFrame(bp_refmir_aln, columns=["bp","refmir"], dtype='Int64')
        self.bp_refmir_aln.set_index("bp", inplace=True)

        return True

    def get_xy(self, i):
        df = self.aln.loc[i]
        return (df['sample'], df['refmir']-0.5)

    def plot_scatter(self, ax, real_start=False, samp_min=None, samp_max=None):
        if samp_min is None: samp_min = 0
        if samp_max is None: samp_max = self.aln['sample'].max()
        i = (self.aln['sample'] >= samp_min) & (self.aln['sample'] <= samp_max)

        return ax.scatter(self.aln['sample'][i], self.aln['refmir'][i], color='orange', zorder=2,s=20)

    def plot_step(self, ax, real_start=False, samp_min=None, samp_max=None):
        i = (self.aln['sample'] >= samp_min) & (self.aln['sample'] <= samp_max)

        ret = ax.step(self.aln['sample'][i], self.aln['refmir'][i], color='orange', zorder=1, where='post')

        if self.errs is not None:
            for t in self.ERR_TYPES:
                e = self.errs[self.errs['type'] == t]
                ax.scatter(
                    e['sample'], e['refmir'], 
                    color='red', zorder=3, 
                    s=self.ERR_SIZES[t],
                    linewidth=self.ERR_WIDTHS[t],
                    marker=self.ERR_MARKS[t]
                )

        return ret

