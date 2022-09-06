#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict, namedtuple
import re
import time
import pandas as pd
import scipy.stats
import copy

from ..pafstats import parse_paf, PafEntry
from ..config import Config
from ..argparse import Opt
from ..index import RefCoord

class Bcaln:
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


    def __init__(self, conf, ref_index, read, paf, clip_coords=None):

        self.is_rna = conf.is_rna

        self.clip_coords = clip_coords

        #ref_coord = RefCoord(paf.rf_name, paf.rf_st-1, paf.rf_en+2, paf.is_fwd)
        #ref_coord = RefCoord(paf.rf_name, paf.rf_st, paf.rf_en, paf.is_fwd)
        self.is_fwd = not paf.is_reverse

        ref_coord = RefCoord(paf.reference_name, paf.reference_start, paf.reference_end, self.is_fwd)
        self.paf_coords = ref_index.get_coord_space(ref_coord, self.is_rna, load_kmers=False)

        self.kmer_shift = ref_index.trim#[not paf.is_fwd]
        #if paf.is_fwd == self.is_rna:
        #self.kmer_shift = self.kmer_shift[::-1]

        self.ref_gaps = list()
        self.errors = None

        self.flip_ref = self.is_fwd == self.is_rna

        if not read.bc_loaded or (not self.parse_cs(paf) and not self.parse_cigar(paf)):
            return

        #TODO make c++ build this 
        moves = np.array(read.moves, bool)
        bce_qrs = np.cumsum(read.moves)
        bce_samps = read.template_start + np.arange(len(bce_qrs)) * read.bce_stride

        samp_bps = pd.DataFrame({
            "start" : bce_samps,
            "length" : read.bce_stride,
            "bp"     : np.cumsum(read.moves),
        })

        df = samp_bps.join(self.bp_mref_aln, on="bp").dropna()

        grp = df.groupby("mref")

        df = pd.DataFrame({
            "mref"    : grp["mref"].first().astype("int64"),
            "start"  : grp["start"].min().astype("uint32"),
            "length" : grp["length"].sum().astype("uint32"),
            "indel" : grp["indel"].first()
        }).set_index("mref")

        #print((df["indel"] < 0).mean(), (df["indel"] > 0).mean())

        #df = pd.concat([df, self.errors], axis=1)

        if self.clip_coords is not None:
            mrefs = df.index.intersection(self.clip_coords.mrefs[self.is_fwd])
            mrefs.name = "mref"

            df = df.reindex(index=mrefs, copy=False)
            self.coords = self.clip_coords#.mref_intersect(mrefs=self.df.index)
        else:
            self.coords = self.paf_coords


        df = df.set_index(df.index - self.kmer_shift[0])

        self.df = df.iloc[self.kmer_shift[0]:]#:-self.kmer_shift[1]]
        self.coords = self.coords.mref_intersect(mrefs=self.df.index)


    @property
    def empty(self):
        return not hasattr(self, "df") or len(self.df) <= sum(self.kmer_shift)

    def parse_cs(self, paf):
        if not paf.has_tag("cs"): return False
        cs = paf.get_tag('cs', (None,)*2)[0]

        #TODO rename to general cig/cs
        bp_mref_aln = list()
        errors = list()

        if not self.is_rna:
            read_i = paf.query_start
        else:
            read_i = paf.query_length - paf.query_end

        mrefs = self.paf_coords.mrefs# - self.kmer_shift[0]
        mref_i = mrefs.min()

        cs_ops = re.findall("(=|:|\*|\+|-|~)([A-Za-z0-9]+)", cs)

        if self.flip_ref:
            cs_ops = reversed(cs_ops)

        for op in cs_ops:
            c = op[0]
            if c in {'=',':'}:
                l = len(op[1]) if c == '=' else int(op[1])
                for qr, mr in zip(range(read_i, read_i+l), range(mref_i, mref_i+l)):
                    if mr in mrefs:
                        bp_mref_aln.append((qr,mr))
                read_i += l
                mref_i += l
            else:
                errors.append( (mref_i,"".join(op)) )

                if c == '*':
                    bp_mref_aln.append((read_i,mref_i))
                    read_i += 1
                    mref_i += 1

                elif c == '-':
                    l = len(op[1])
                    mref_i += l

                elif c == '+':
                    l = len(op[1])
                    read_i += l

                elif c == '~':
                    l = int(op[1][2:-2])
                    self.ref_gaps.append((mref_i,mref_i+l))
                    mref_i += l

                else:
                    print("UNIMPLEMENTED ", op)

        self.bp_mref_aln = pd.DataFrame(bp_mref_aln, columns=["bp","mref"], dtype='Int64')
        self.bp_mref_aln.set_index("bp", inplace=True)

        #TODO type shouldn't have to be 64 bit
        self.errors = pd.DataFrame(errors, columns=["mref","error"]) \
                       .set_index("mref").groupby(level=0) \
                       .transform(lambda errs: ",".join(errs))

        return True        

    def parse_cigar(self, paf):
        #cig = paf.tags.get('cg', (None,)*2)[0]
        cig = paf.cigarstring
        if cig is None: return False

        bp_mref_aln = list()#defaultdict(list)

        #print(paf.query_alignment_start, paf.query_alignment_end)
        #if not self.is_rna:
        #    read_i = paf.query_alignment_start
        #else:
        #    read_i = paf.infer_query_length() - paf.query_alignment_end
        read_i = 0

        cig_ops = self.CIG_RE.findall(cig)

        if self.is_fwd == self.is_rna:
            cig_ops = list(reversed(cig_ops))

        mrefs = self.paf_coords.mrefs# - self.kmer_shift[0]
        mref_i = mrefs.min()

        insert_len = 0

        for l,c in cig_ops:
            l = int(l)
            incr_qr = c in self.CIG_INCR_RD
            incr_rf = c in self.CIG_INCR_RF
            read_j = read_i + (l if incr_qr else 1)
            mref_j = mref_i + (l if incr_rf else 1)

            #CIG_INCR_RD = CIG_INCR_ALL | {'I','S'}
            #CIG_INCR_RF = CIG_INCR_ALL | {'D','N'}

            if c == "M":
                for qr, mr in zip(range(read_i, read_j), range(mref_i, mref_j)):
                    #if mr in mrefs:
                    bp_mref_aln.append((qr,mr,insert_len))
                    insert_len = 0
            elif c == "I":
                insert_len = l

            elif c in {'D','N'}:
                for rf in range(mref_i, mref_j):
                    bp_mref_aln.append((read_i,rf,-l))
                    
                #TODO do this based on indel field
                if c == "N":
                    self.ref_gaps.append((mref_i,mref_j))
                
            if incr_qr: read_i = read_j 
            if incr_rf: mref_i = mref_j 

        self.bp_mref_aln = pd.DataFrame(bp_mref_aln, columns=["bp","mref","indel"], dtype='Int64')
        self.bp_mref_aln.set_index("bp", inplace=True)

        return True

