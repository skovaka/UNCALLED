#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import division
import sys                         
import os
import numpy as np
from bisect import bisect_left, bisect_right
from typing import NamedTuple
import collections.abc
import pysam

import pandas as pd
from . import RefCoord, str_to_coord, SeqRecord

import _uncalled
from .pore_model import PoreModel
from .argparse import Opt

class FastaIndex:
    def __init__(self, model, filename):
        self.infile = pysam.FastaFile(filename)
        self.model = model

        if isinstance(self.model.instance, _uncalled.PoreModelU16):
            self.SeqInst = _uncalled.SequenceU16
        elif isinstance(self.model.instance, _uncalled.PoreModelU32):
            self.SeqInst = _uncalled.SequenceU32
        else:
            raise ValueError(f"Unknown PoreModel type: {model.instance}")

        self.ref_ids = dict()
        self.anno = list()
        self.offsets = list()
        offs = 0
        for i in range(self.infile.nreferences):
            name = self.infile.references[i]
            size = self.infile.lengths[i]
            self.anno.append(SeqRecord(name, i, size, offs))
            self.offsets.append(offs)
            offs += size

    def pac_to_pos(self, pac):
        if hasattr(pac, "__len__"):
            c = pac[0]
        else:
            c = pac
        i = np.searchsorted(self.offsets, c, side="right")-1
        if i not in range(len(self.offsets)):
            raise ValueError(f"Packed sequence coordinate out of range: {c}")
        return pac - self.offsets[i]

    def pos_to_mpos(self, pos, fwd):
        if fwd == self.model.reverse:
            return -pos-1
        return pos

    def mpos_to_pos(self, mpos):
        if hasattr("__len__", pac):
            c = pac[0]
        else:
            c = pac
        if c < 0:
            return -mpos-1;
        return mpos

    def query(self, coord):
        seqs = list()
        bounds = coord.bounds
        for i in range(0,len(bounds),2):
            seqs.append(self.infile.fetch(coord.name, bounds[i], bounds[i+1]))
        return self.SeqInst(self.model.instance, "".join(seqs), coord)

    def get_ref_name(self, rid):
        return self.anno[rid].name

    def get_pac_offset(self, rid):
        return self.anno[rid].offset

_index_cache = dict()

def load_index(model, prefix, load_pacseq=True, load_bwt=False, cache=True):
    idx = _index_cache.get(prefix, None)
    if idx is None:
#        idx = BwaIndex(model, prefix, load_pacseq, load_bwt)
        idx = FastaIndex(model, prefix)
        if cache: _index_cache[prefix] = idx
    return idx

