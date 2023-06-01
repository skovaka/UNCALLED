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
from ..pore_model import PoreModel
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
        return None
        
    return moves_to_aln(moves, template_start, mv_stride)

def sam_to_ref_moves(conf, ref_index, read, sam, workflow=None):
    if read is None or read.empty(): 
        return None
    read_moves = sam_to_read_moves(read, sam)
    if read_moves is None:
        return None

    model = ref_index.model

    is_fwd = not sam.is_reverse
    flip_ref = is_fwd == model.reverse

    ar = np.array(sam.get_aligned_pairs())
    ar = ar[ar[:,1] != None] #TODO keep track of insertion counts


    if flip_ref:
        ar = ar[::-1]
        qrys = ar[:,0]
        #silly trick to make null reverse into REF_NA
        read_len = sam.infer_read_length()
        qrys[qrys == None] = read_len - INT32_NA - 1
        qrys = read_len - qrys - 1
    else:
        qrys = ar[:,0]
        qrys[qrys == None] = INT32_NA

    refs = np.array(ref_index.pos_to_mpos(ar[:,1], is_fwd))
    qrys = np.array(qrys, dtype=np.int64)


    ref_moves = read_to_ref_moves(read_moves, refs, qrys, conf.dtw.del_max, True)

    if workflow is None:
        workflow = (conf.read_buffer.flowcell, conf.read_buffer.kit)

    mkl = MOVE_KMER_LENS[PoreModel.PRESET_MAP.loc[workflow, "ont_model"]]

    shift = model.K - mkl - model.shift
    ref_moves.index.shift(shift)

    return ref_moves.slice(-shift+model.shift, len(ref_moves)-model.K+model.shift+1)
