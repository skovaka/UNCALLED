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

from ..config import Config
from ..argparse import Opt
from ..index import RefCoord

class Bcaln:
    CIG_OPS_STR = "MIDNSHP=X"
    CIG_RE = re.compile("(\d+)(["+CIG_OPS_STR+"])")
    CIG_OPS = set(CIG_OPS_STR)
    CIG_INCR_ALL = {'M','=', 'X'}
    CIG_INCR_RD = CIG_INCR_ALL | {'I','S','H'}
    CIG_INCR_RF = CIG_INCR_ALL | {'D','N'}

    SUB = 0
    INS = 1
    DEL = 2
    ERR_TYPES = [SUB, INS, DEL]
    ERR_MARKS = ['o', 'P', '_']
    ERR_SIZES = [100, 150, 150]
    ERR_WIDTHS = [0,0,5]


    def __init__(self, conf, ref_index, read, sam, clip_coords=None):

        self.is_rna = conf.is_rna
        self.del_max = conf.dtw.del_max
        self.ins_max = conf.dtw.ins_max

        self.clip_coords = clip_coords

        #ref_coord = RefCoord(sam.rf_name, sam.rf_st-1, sam.rf_en+2, sam.is_fwd)
        #ref_coord = RefCoord(sam.rf_name, sam.rf_st, sam.rf_en, sam.is_fwd)
        self.is_fwd = not sam.is_reverse

        ref_coord = RefCoord(sam.reference_name, sam.reference_start, sam.reference_end, self.is_fwd)
        self.sam_coords = ref_index.get_coord_space(ref_coord, self.is_rna, load_kmers=False)

        self.flip_ref = self.is_fwd == self.is_rna

        self.kmer_shift = ref_index.trim#[not sam.is_fwd]
        if not self.flip_ref:
            self.kmer_shift = self.kmer_shift[::-1]
        #if sam.is_fwd == self.is_rna:
        #self.kmer_shift = self.kmer_shift[::-1]

        self.ref_gaps = list()
        self.errors = None

        if read.bc_loaded:
            moves = read.moves
        elif sam.has_tag("mv"):
            moves = np.array(sam.get_tag("mv"))[1:]
        else:
            raise RuntimeError(f"Basecaller \"moves\" not found for read {read.id}")
            
        if not self.parse_cigar(sam):
            raise RuntimeError(f"Cigar string not found for read {read.id}")

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
        #pd.set_option('display.max_rows', 50000) 
        #print(df)

        #print((df["indel"] < 0).mean(), (df["indel"] > 0).mean())

        #df = pd.concat([df, self.errors], axis=1)

        if self.clip_coords is not None:
            mrefs = df.index.intersection(self.clip_coords.mrefs[self.is_fwd])
            mrefs.name = "mref"

            df = df.reindex(index=mrefs, copy=False)
            self.coords = self.clip_coords#.mref_intersect(mrefs=self.df.index)
        else:
            self.coords = self.sam_coords

        #TODO don't do this
        if self.kmer_shift[0] <= 2:
            df = df.set_index(df.index - self.kmer_shift[0])
        else:
            df = df.set_index(df.index - self.kmer_shift[0] + 1)

        self.df = df.iloc[self.kmer_shift[0]:]#:-self.kmer_shift[1]]
        self.coords = self.coords.mref_intersect(mrefs=self.df.index)


    @property
    def empty(self):
        return not hasattr(self, "df") or len(self.df) <= sum(self.kmer_shift)

    def parse_cigar(self, sam):
        #cig = sam.tags.get('cg', (None,)*2)[0]
        cig = sam.cigarstring
        if cig is None: return False

        bp_mref_aln = list()#defaultdict(list)

        #print(sam.query_alignment_start, sam.query_alignment_end)
        #if not self.is_rna:
        #    read_i = sam.query_alignment_start
        #else:
        #    read_i = sam.infer_query_length() - sam.query_alignment_end
        read_i = 0

        cig_ops = self.CIG_RE.findall(cig)

        if self.is_fwd == self.is_rna:
            cig_ops = list(reversed(cig_ops))

        mrefs = self.sam_coords.mrefs# - self.kmer_shift[0]
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
                if c == "N" or l > self.del_max:
                    self.ref_gaps.append((mref_i,mref_j))
                
            if incr_qr: read_i = read_j 
            if incr_rf: mref_i = mref_j 

        self.bp_mref_aln = pd.DataFrame(bp_mref_aln, columns=["bp","mref","indel"], dtype='Int64')
        self.bp_mref_aln.set_index("bp", inplace=True)

        return True

