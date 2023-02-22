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

from _uncalled import _AlnDF, IntervalIndexI64, IntervalIndexI32, moves_to_aln, read_to_ref_moves
from ..config import Config
from ..argparse import Opt
from ..index import RefCoord

MOVE_KMER_LENS = {
   "dna_r10.3_450bps" : 10,
   "dna_r10.3_450bpsm" : 10,
   "dna_r10.4.1_e8.2_260bps" : 10,
   "dna_r10.4.1_e8.2_400bps" : 10,
   "dna_r10_450bps" : 10,
   "dna_r10.4_e8.1" : 10,
   "dna_r10.4_e8.1m" : 10,
   "dna_r9.4.1_450bps" : 6,
   "dna_r9.4.1_e8.1" : 6,
   "dna_r9.4.1_e8.1m" : 6,
   "dna_r9.5_450bps" : 6,
   "rna_r9.4.1_70bps" : 5,
}

class Bcaln:
    CIG_OPS_STR = "MIDNSHP=X"
    CIG_RE = re.compile("(\d+)(["+CIG_OPS_STR+"])")
    CIG_OPS = set(CIG_OPS_STR)
    CIG_INCR_ALL = {'M','=', 'X'}
    CIG_INCR_RD = CIG_INCR_ALL | {'I','S', 'H'}
    CIG_INCR_RF = CIG_INCR_ALL | {'D','N'}

    SUB = 0
    INS = 1
    DEL = 2
    ERR_TYPES = [SUB, INS, DEL]
    ERR_MARKS = ['o', 'P', '_']
    ERR_SIZES = [100, 150, 150]
    ERR_WIDTHS = [0,0,5]


    def __init__(self, conf, ref_index, read, sam, clip_coords=None, model_name=None):

        self.is_rna = conf.is_rna
        self.del_max = conf.dtw.del_max
        self.ins_max = conf.dtw.ins_max

        self.clip_coords = clip_coords

        if model_name is not None:
            mkl = MOVE_KMER_LENS[model_name]
        else:
            mkl = conf.pore_model.k

        self.is_fwd = not sam.is_reverse

        ref_coord = RefCoord(sam.reference_name, sam.reference_start, sam.reference_end, self.is_fwd)
        self.sam_coords = ref_index.get_coord_space(ref_coord, self.is_rna, load_kmers=False)

        self.flip_ref = self.is_fwd == self.is_rna

        self.kmer_shift = ref_index.trim#[not sam.is_fwd]
        if self.flip_ref:
            self.kmer_shift = self.kmer_shift[::-1]

        self.ref_gaps = list()
        self.errors = None

        if read.bc_loaded:
            mv_stride = read.bce_stride
            moves = np.array(read.moves)
            template_start = read.template_start

        elif sam.has_tag("mv"):
            mv = np.array(sam.get_tag("mv"))
            mv_stride = mv[0]
            template_start = sam.get_tag("ts")
            moves = mv[1:]

        else:
            sys.stderr.write(f"Basecaller moves not found for read {read.id}, skipping\n")
            self.df = pd.DataFrame()
            return
            
        if not self.parse_cigar(sam):
            raise RuntimeError(f"Cigar string not found for read {read.id}")

        #df = pd.DataFrame(
        #    sam.get_aligned_pairs(), 
        #    columns=["read","ref"], dtype="Int64"
        #).dropna(subset=["ref"])

        #df["mref"] = self.sam_coords.ref_to_mref(df["ref"].to_numpy())

        #if self.flip_ref:
        #    df["read"] = sam.query_length-df["read"]-1
        #    df = df.iloc[::-1]


        #self.bp_mref_aln = df.set_index("read")["mref"]
        # .set_index("ref")\
        # .reindex(self.sam_coords.refs)\
        # .reset_index()\
        # .set

        #df = df.set_index(self.sam_coords.ref_to_mref(df.index))
        #print(self.bp_mref_aln)

        #TODO make c++ build this 
        #moves = np.array(read.moves, bool)
        bce_qrs = np.cumsum(moves)
        bce_samps = template_start + np.arange(len(bce_qrs)) * mv_stride

        read_moves = moves_to_aln(moves, template_start, mv_stride)

        QRY_NA = np.iinfo(np.int32).max

        ar = np.array(sam.get_aligned_pairs())
        ar = ar[ar[:,1] != None]
        if self.flip_ref:
            ar = ar[::-1]
            qrys = ar[:,0]
            #silly trick to make null reverse into REF_NA
            qrys[qrys == None] = sam.query_length - QRY_NA - 1
            qrys = sam.query_length - qrys - 1

        else:
            qrys = ar[:,0]
            qrys[qrys == None] = QRY_NA
        refs = np.array(self.sam_coords.ref_to_mref(ar[:,1]), np.int64)
        qrys = qrys.astype(np.int64)

        ref_moves = read_to_ref_moves(read_moves, refs, qrys, conf.dtw.del_max)
        df= pd.DataFrame({
            "mref"   : np.array(ref_moves.index.expand()),
            "start"  : np.array(ref_moves.samples.starts),
            "length" : np.array(ref_moves.samples.lengths)
        }).set_index("mref")

        isna = df["start"] == QRY_NA
        df[isna] = pd.NA
        df.fillna(method="backfill", inplace=True)

        #ref_moves = RefMoves(read_moves, refs, qrys)

        #samp_bps = pd.DataFrame({
        #    "start" : bce_samps,
        #    "length" : mv_stride,
        #    "bp"     : bce_qrs,
        #})

        #samp_bps = samp_bps.join(self.bp_mref_aln, on="bp").dropna()
        #grp = samp_bps.groupby("mref")

        #df = pd.DataFrame({
        #    "mref"    : grp["mref"].first().astype("int64"),
        #    "start"  : grp["start"].min().astype("uint32"),
        #    "length" : grp["length"].sum().astype("uint32"),
        #    #"indel" : grp["indel"].first()
        #}).set_index("mref")

        #print("sorted")
        #print(df.sort_values("length"))

        #print("diff")
        #print(df[diff])
        #print()

        if self.clip_coords is not None:
            mrefs = df.index.intersection(self.clip_coords.mrefs[self.is_fwd])
            mrefs.name = "mref"

            df = df.reindex(index=mrefs, copy=False)
            self.coords = self.clip_coords#.mref_intersect(mrefs=self.df.index)
        else:
            self.coords = self.sam_coords

        #TODO don't do this
        #if self.kmer_shift[0] <= 2:
        #if self.kmer_shift[0] <= 2:
        #    df = df.set_index(df.index - self.kmer_shift[0])
        #else:
        #    df = df.set_index(df.index - self.kmer_shift[0] + 1)

        sh = conf.pore_model.k - mkl
        self.df = df.set_index(df.index - (self.kmer_shift[0] - sh))\
                    .reindex(self.coords.mrefs[self.kmer_shift[0]:-self.kmer_shift[1]].sort_values() + sh)
        self.df.index.name = "mref"
        #df = df.set_index(df.index - (self.kmer_shift[0] - sh))
        #self.df = df.iloc[np.sum(self.kmer_shift) - sh:len(df)-sh]#-self.kmer_shift[1]+sh]
        self.coords = self.coords.mref_intersect(mrefs=self.df.index)


    @property
    def empty(self):
        return not hasattr(self, "df") or len(self.df) <= sum(self.kmer_shift)

    def parse_cigar(self, sam):
        #cig = sam.tags.get('cg', (None,)*2)[0]
        cig = sam.cigarstring
        if cig is None: return False

        bp_mref_aln = list()#defaultdict(list)

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

