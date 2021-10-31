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
from ..config import Config
from ..argparse import Opt
from .. import nt, PoreModel
from ..index import RefCoord

class Bcaln:
    K = 4
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


    def __init__(self, ref_index, read, paf, clip_coords=None):

        self.is_rna = read.conf.is_rna

        self.clip_coords = clip_coords

        ref_coord = RefCoord(paf.rf_name, paf.rf_st, paf.rf_en, paf.is_fwd)
        self.paf_coords = ref_index.get_coord_space(ref_coord, self.is_rna, kmer_trim=True)

        self.refgap_bps = list()
        self.sub_bps = list()
        self.ins_bps = list()
        self.del_bps = list()
        self.err_bps = None

        self.is_fwd = paf.is_fwd
        self.flip_ref = paf.is_fwd == self.is_rna

        if not read.f5.bc_loaded or (not self.parse_cs(paf) and not self.parse_cigar(paf)):
            return
        #if self.empty: 
        #    return

        #TODO make c++ build this 
        moves = np.array(read.f5.moves, bool)
        bce_qrs = np.cumsum(read.f5.moves)
        bce_samps = read.f5.template_start + np.arange(len(bce_qrs)) * read.f5.bce_stride

        samp_bps = pd.DataFrame({
            "start" : bce_samps,
            "length" : read.f5.bce_stride,
            "bp"     : np.cumsum(read.f5.moves),
        })

        df = samp_bps.join(self.bp_mref_aln, on="bp").dropna()
        #df["mref"] = df["mref"].astype("Int64")

        grp = df.groupby("mref")

        df = pd.DataFrame({
            "mref"    : grp["mref"].first().astype("int64"),
            "start"  : grp["start"].min().astype("uint32"),
            "length" : grp["length"].sum().astype("uint32"),
            "bp"   : grp["bp"].first()
        }).set_index("mref")

        if self.err_bps is not None:
            self.errs = samp_bps.join(self.err_bps.set_index("bp"), on="bp").dropna()
            self.errs.reset_index(inplace=True, drop=True)
        else:
            self.errs = None

        if self.clip_coords is not None:
            mrefs = df.index.intersection(self.clip_coords.mrefs[self.is_fwd])

            self.coords = self.clip_coords.mref_intersect(mrefs=df.index)

            #mrefs = df.index.intersect(self.clip_coords.mrefs[self.is_fwd])
            df = df.reindex(index=mrefs, copy=False)
        else:
            self.coords = self.paf_coords

        self.df = df


    @property
    def empty(self):
        return not hasattr(self, "df") or len(self.df) <= 2

    def parse_cs(self, paf):
        cs = paf.tags.get('cs', (None,)*2)[0]
        if cs is None: return False

        #TODO rename to general cig/cs
        bp_mref_aln = list()
        err_bps = list()

        if not self.is_rna:
            qr_i = paf.qr_st
            #rf_i = paf.rf_st
        else:
            qr_i = paf.qr_len - paf.qr_en 
            #rf_i = -paf.rf_en+1

        mrefs = self.paf_coords.mrefs
        mr_i = mrefs.min()

        cs_ops = re.findall("(=|:|\*|\+|-|~)([A-Za-z0-9]+)", cs)

        if self.flip_ref:
            cs_ops = reversed(cs_ops)

        for op in cs_ops:
            c = op[0]
            if c in {'=',':'}:
                l = len(op[1]) if c == '=' else int(op[1])
                for qr, mr in zip(range(qr_i, qr_i+l), range(mr_i, mr_i+l)):
                    if mr in mrefs:
                        bp_mref_aln.append((qr,mr))
                #bp_mref_aln += zip(range(qr_i, qr_i+l), range(mr_i, mr_i+l))
                qr_i += l
                mr_i += l

            elif c == '*':
                self.sub_bps.append(qr_i)
                bp_mref_aln.append((qr_i,mr_i))
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

        self.bp_mref_aln = pd.DataFrame(bp_mref_aln, columns=["bp","mref"], dtype='Int64')
        self.bp_mref_aln.set_index("bp", inplace=True)

        #TODO type shouldn't have to be 64 bit
        self.err_bps = pd.DataFrame(err_bps, columns=["bp","mref","type","seq"])#, dtype='Int64')

        return True        

    def parse_cigar(self, paf):
        cig = paf.tags.get('cg', (None,)*2)[0]
        if cig is None: return False

        bp_mref_aln = list()#defaultdict(list)
        self.refgap_bps = list()

        #mr_i = self.mref_start
        if not self.is_rna:
            qr_i = paf.qr_st
        else:
            qr_i = paf.qr_len - paf.qr_en 

        mrefs = self.paf_coords.mrefs
        mr_i = mrefs.min()

        cig_ops = self.CIG_RE.findall(cig)

        if paf.is_fwd == self.is_rna:
            cig_ops = list(reversed(cig_ops))

        for l,c in cig_ops:
            l = int(l)
            incr_qr = c in self.CIG_INCR_RD
            incr_rf = c in self.CIG_INCR_RF
            qr_j = qr_i + (l if incr_qr else 1)
            mr_j = mr_i + (l if incr_rf else 1)

            if c == "M":
                for qr, mr in zip(range(qr_i, qr_j), range(mr_i, mr_j)):
                    if mr in mrefs:
                        bp_mref_aln.append((qr,mr))
            elif c == "N":
                if mr_i in mrefs:
                    bp_mref_aln.append((qr_i,mr))

            if incr_qr:
                qr_i = qr_j 

            if incr_rf:
                mr_i = mr_j 

        self.bp_mref_aln = pd.DataFrame(bp_mref_aln, columns=["bp","mref"], dtype='Int64')
        self.bp_mref_aln.set_index("bp", inplace=True)

        return True

