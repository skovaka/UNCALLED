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
INT32_NA = np.iinfo(np.int32).max

def sam_to_read_moves(read, sam):
    if read.bc_loaded:
        mv_stride = read.move_stride
        moves = np.array(read.moves)
        template_start = read.template_start

    elif sam.has_tag("mv"):
        mv = np.array(sam.get_tag("mv"))
        mv_stride = mv[0]
        template_start = sam.get_tag("ts")
        moves = mv[1:]

    else:
        sys.stderr.write(f"Basecaller moves not found for read {read.id}, skipping\n")
        return None
        
    return moves_to_aln(moves, template_start, mv_stride)

def sam_to_ref_moves(conf, ref_index, read, sam, model_name=None):
    if read is None: 
        return None
    read_moves = sam_to_read_moves(read, sam)
    if read_moves is None:
        return None

    is_fwd = not sam.is_reverse
    flip_ref = is_fwd == conf.is_rna

    ar = np.array(sam.get_aligned_pairs())
    ar = ar[ar[:,1] != None] #TODO keep track of insertion counts

    if flip_ref:
        ar = ar[::-1]
        qrys = ar[:,0]
        #silly trick to make null reverse into REF_NA
        qrys[qrys == None] = sam.query_length - INT32_NA - 1
        qrys = sam.query_length - qrys - 1
    else:
        qrys = ar[:,0]
        qrys[qrys == None] = INT32_NA

    refs = ref_index.pos_to_mpos(sam.reference_id, ar[:,1], is_fwd, conf.is_rna)
    qrys = qrys.astype(np.int64)

    ref_moves = read_to_ref_moves(read_moves, refs, qrys, conf.dtw.del_max, True)

    if model_name is not None:
        mkl = MOVE_KMER_LENS[model_name]
    else:
        mkl = conf.pore_model.k

    shift = ref_index.K - mkl - ref_index.trim[0]
    ref_moves.index.shift(shift)
    
    return ref_moves.slice(-shift+ref_index.trim[0], len(ref_moves)-ref_index.trim[1])


class Bcaln:

    def __init__(self, conf, ref_index, read, sam, clip_coords=None, model_name=None):
        self.aln = sam_to_ref_moves(conf, ref_index, read, sam, model_name)

        if self.aln is None:
            self.df = pd.DataFrame()
            return

        df = pd.DataFrame({
            "mpos"   : (self.aln.index.expand()),
            "start"  : (self.aln.samples.starts),
            "length" : (self.aln.samples.lengths)
        }).set_index("mpos")
        isna = df["start"] == INT32_NA
        df[isna] = pd.NA
        self.df = df.fillna(method="backfill")

        self.df.index.name = "mpos"

        ref_coord = RefCoord(sam.reference_name, sam.reference_start, sam.reference_end, not sam.is_reverse)
        sam_coords = ref_index.get_coord_space(ref_coord, conf.is_rna, load_kmers=False)
        self.coords = sam_coords.mpos_intersect(mposs=self.df.index)
        self.K  = conf.pore_model.k

    @property
    def empty(self):
        return not hasattr(self, "df") or len(self.df) <= self.K
