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
import uncalled as unc
from bisect import bisect_left, bisect_right
from typing import NamedTuple
import collections.abc

import pandas as pd

from _uncalled import _RefIndex
from .argparse import Opt
from .config import ParamGroup
from . import nt

class RefIndex(_RefIndex):
    class CoordSpace:
        def __init__(self, ref_name, refs, fwd=None, mrefs=None, index=None, is_rna=None, load_kmers=False):
            self.ref_name = ref_name

            self.refs = refs#pd.RangeIndex(ref_bounds.start, ref_bounds.end-nt.K+1)
            self.fwd = fwd

            if mrefs is None:
                if self.fwd is None:
                    self.mrefs = tuple( (self._init_mrefs(index, fwd, is_rna) for fwd in range(2)) )
                else:
                    self.mrefs = self._init_mrefs(index, self.fwd, is_rna)
            else:
                self.mrefs = mrefs

            if isinstance(self.mrefs, pd.Index):
                self.stranded = True
            elif isinstance(self.mrefs, tuple):
                self.stranded = False
            else:
                raise ValueError("Coord space mref can only be initialized using pandas.RangeIndex or tuple of two RangeIndexes")

            if load_kmers:
                self.kmers = self.load_kmers(index)
            else:
                self.kmers = None
            
        def _init_mrefs(self, index, fwd, is_rna):
            st,en = index.ref_to_mref(self.ref_name, self.refs.start, self.refs.stop, fwd, is_rna)
            mrefs = pd.RangeIndex(st,en)
            if index.is_mref_flipped(mrefs.start):
                return mrefs[::-1]
            return mrefs

        def load_kmers(self, index):
            if self.stranded:
                return self._mrefs_to_kmers(index, self.mrefs, self.fwd)
            return tuple( (self._mrefs_to_kmers(index, mrefs, fwd) 
                            for fwd,mrefs in enumerate(self.mrefs)) )

        def _mrefs_to_kmers(self, index, mrefs, fwd):
            kmers = index.get_kmers(mrefs.min(), mrefs.max()+nt.K, fwd)
            if self.mrefs.step < 0:
                kmers = kmers[::-1]
            return pd.Series(index=mrefs, data=kmers)

        def intersect(self, coords):
            if self.ref_name != coords.ref_name:
                raise ValueError("Cannot intersect CoordSpaces from different reference sequences")

            ref_min = max(self.refs.min(), coords.refs.min())
            ref_max = min(self.refs.max(), coords.refs.max())
            ref_len = ref_max - ref_min + 1

            st = self.refs.get_loc(ref_min)
            en = st + ref_len

            refs = self.refs[st:en]#self.refs.intersection(coords.refs)

            if self.stranded:
                mrefs = self.mrefs[st:en] #.intersection(coords.mrefs)
            else:
                mrefs = tuple( (m[st:en] for m in self.mrefs) )

            print(refs, mrefs)
            return RefIndex.CoordSpace(self.ref_name, refs, mrefs=mrefs)

        def ref_to_mref(self, ref, fwd=None):
            if isinstance(ref, (collections.abc.Sequence, np.ndarray)):
                i = self.refs.get_indexer(ref)
            else:
                i = self.refs.get_loc(ref)

            if self.stranded:
                if fwd is not None and fwd != self.fwd:
                    raise ValueError("Invalid 'fwd' value or stranded CoordSpace")
                return self.mrefs[i]

            # not stranded
            if fwd is None:
                return tuple( (mrefs[i] for mrefs in self.mrefs) )
            return self.mrefs[fwd][i]

        def is_mref_fwd(self, mref):
            if isinstance(mref, (collections.abc.Sequence, np.ndarray)):
                mref = mref[0]

            if self.stranded:
                return self.fwd and mref in self.mrefs
            else:
                return mref in self.mrefs[1]
                
        #TODO make single private method with fwd param. Also probably merge with above
        def all_mrefs_fwd(self, mrefs):
            return len(self.mrefs[True].intersection(mrefs)) == len(mrefs)
                
        def all_mrefs_rev(self, mrefs):
            return len(self.mrefs[False].intersection(mrefs)) == len(mrefs)
        
        def mref_to_ref(self, mref):
            if self.stranded:
                mrefs = self.mrefs
            else:
                mrefs = self.mrefs[self.is_mref_fwd(mref)]

            if isinstance(mref, (collections.abc.Sequence, np.ndarray)):
                i = mrefs.get_indexer(mref)
            else:
                i = mrefs.get_loc(mref)

            return self.refs[i]

        #def validate_refs(self, refs):
        #    return len(self.refs.intersection(refs)) == len(refs)

        #def validate_mrefs(self, mrefs):
        #    return all_mrefs_fwd(mrefs) or all_mrefs_rev(mrefs)

        def __len__(self):
            return len(self.refs)

        def __repr__(self):
            return "%s:%d-%d (fwd %d-%d, rev %d-%d)" % (
                self.ref_name, self.refs.start, self.refs.stop,
                self.mrefs[True].start, self.mrefs[True].stop,
                self.mrefs[False].start, self.mrefs[False].stop,
            )

    def get_coord_space(self, ref_bounds, is_rna, kmer_shift=nt.K-1):
        rid = self.get_ref_id(ref_bounds.name)
        length = self.get_ref_len(rid)
        if ref_bounds.start < 0 or ref_bounds.end > length:
            raise ValueError("Reference coordinates %s out of bounds for sequence of length %d" % (ref_bounds, length))

        #TODO get of -nt.K+1
        refs = pd.RangeIndex(ref_bounds.start, ref_bounds.end-kmer_shift)
        return self.CoordSpace(ref_bounds.name, refs, fwd=ref_bounds.fwd, index=self, is_rna=is_rna)

_index_cache = dict()

def load_index(prefix, load_pacseq=True, load_bwt=False, cache=True):
    idx = _index_cache.get(prefix, None)
    if idx is None:
        idx = RefIndex(prefix, load_pacseq, load_bwt)
        if cache: _index_cache[prefix] = idx
    else:
        if load_pacseq and not idx.pacseq_loaded():
            idx.load_pacseq()
        if load_bwt and not idx.bwt_loaded():
            idx.load_index()
    return idx

#Index parameter group
class IndexParams(ParamGroup):
    _name = "index"
IndexParams._def_params(
    ("fasta_filename", None, str, "FASTA file to index"),
    ("index_prefix", None, str, "Index output prefix. Will use input fasta filename by default"),
    ("no_bwt", False, bool, "Will only generate the pacseq if specified, which is much faster to build. Can only be used with DTW subcommands (NOT map, sim, or realtime)"),
    ("max_sample_dist", 100, int, "Maximum average sampling distance between reference alignments."),
    ("min_samples", 50000, int, "Minimum number of alignments to produce (approximate, due to deterministically random start locations),"),
    ("max_samples", 1000000, int, "Maximum number of alignments to produce (approximate, due to deterministically random start locations),"),
    ("kmer_len", 5, int, "Model k-mer length"),
    ("matchpr1", 0.6334, float, "Minimum event match probability"),
    ("matchpr2", 0.9838, float, "Maximum event match probability"),
    ("pathlen_percentile", 0.05, float, ""),
    ("max_replen", 100, int, ""),
    ("probs", None, str, "Find parameters with specified target probabilites (comma separated)"),
    ("speeds", None, str, "Find parameters with specified speed coefficents (comma separated)"),
)

BWA_OPTS = (
    Opt("bwa_prefix", "mapper"),
    Opt(("-p", "--idx-preset"), "mapper"),
)

OPTS = (
    Opt("fasta_filename", "index"),
    Opt(("-o", "--index-prefix"), "index"),
    Opt("--no-bwt", "index", action="store_true"),
    Opt(("-s", "--max-sample-dist"), "index"),
    Opt("--min-samples", "index"),
    Opt("--max-samples", "index"),
    Opt(("-1", "--matchpr1"), "index"),
    Opt(("-2", "--matchpr2"), "index"),
    Opt(("-f", "--pathlen-percentile"), "index"),
    Opt(("-m", "--max-replen"), "index"),
    Opt("--probs", "index"),
    Opt("--speeds", "index"),
)

def main(conf):
    """Build an index from a FASTA reference"""
        
#       Uses BWA to generate pacseq and FM index files, then computes reference-specific probability thresholds. The FM index and probability thresholds are only required for Real-Time Enrichment subcommands (see --no-bwt if only interested in DTW analysis)"""

    prms = conf.index

    if prms.index_prefix is None or len(prms.index_prefix) == 0:
        prms.index_prefix = prms.fasta_filename

    bwa_built = True

    for suff in unc.index.UNCL_SUFFS:
        if not os.path.exists(prms.index_prefix + suff):
            bwa_built = False
            break

    if bwa_built:
        sys.stderr.write("Using previously built BWA index.\nNote: to fully re-build the index delete files with the \"%s.*\" prefix.\n" % prms.index_prefix)
    else:
        unc.RefIndex.create(prms.fasta_filename, prms.index_prefix, prms.no_bwt)
        
        if prms.no_bwt: 
            sys.stderr.write("Pacseq built\n")
            return

    sys.stderr.write("Initializing parameter search\n")
    p = unc.index.IndexParameterizer(prms)

    p.add_preset("default", tgt_speed=115)

    if prms.probs != None:
        for tgt in prms.probs.split(","):
            sys.stderr.write("Writing 'prob_%s' parameters\n" % tgt)
            try:
                p.add_preset("prob_%s" % tgt, tgt_prob=float(tgt))
            except Exception as e:
                sys.stderr.write("Failed to add 'prob_%s'\n" % tgt)

    if prms.speeds != None:
        for tgt in prms.speeds.split(","):
            sys.stderr.write("Writing 'speed_%s' parameters\n" % tgt)
            try:
                p.add_preset("speed_%s" % tgt, tgt_speed=float(tgt))
            except:
                sys.stderr.write("Failed to add 'speed_%s'\n" % tgt)

    p.write()

    sys.stderr.write("Done\n")

#TODO move to RefIndex?
UNCL_SUFF = ".uncl"
AMB_SUFF = ".amb"
ANN_SUFF = ".ann"
BWT_SUFF = ".bwt"
PAC_SUFF = ".pac" 
SA_SUFF = ".sa"
NOBWT_SUFFS = [ANN_SUFF, AMB_SUFF, PAC_SUFF]
UNCL_SUFFS = NOBWT_SUFFS + [UNCL_SUFF, PAC_SUFF, SA_SUFF]

def check_prefix(path, no_bwt=False):
    if no_bwt:
        suffs = NOBWT_SUFFS
    else:
        suffs = UNCL_SUFFS

    for suff in suffs:
        fname = path+suff
        if not os.path.exists(fname):
            raise FileNotFoundError("could not find index file \"%s\"" % fname)
    return True

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


def power_fn(xmax, ymin, ymax, exp, N=100):
    dt = 1.0/N
    t = np.arange(0, 1+dt, dt)

    return t*xmax, (t**exp) * (ymax-ymin) + ymin

class IndexParameterizer:
    #MODEL_THRESHS_FNAME = os.path.join(ROOT_DIR, "conf/r94_5mers_rna_threshs.txt")
    MODEL_THRESHS_FNAME = os.path.join(ROOT_DIR, "config/r94_5mers_threshs.txt")

    def __init__(self, params):
        self.prms = params

        self.out_fname = self.prms.index_prefix + UNCL_SUFF

        self.pck1 = self.prms.matchpr1
        self.pck2 = self.prms.matchpr2

        self.calc_map_stats()
        self.get_model_threshs()

        self.functions = dict()

    def calc_map_stats(self):

        ann_in = open(self.prms.index_prefix + ANN_SUFF)
        header = ann_in.readline()
        ref_len = int(header.split()[0])
        ann_in.close()

        approx_samps = ref_len / self.prms.max_sample_dist
        if approx_samps < self.prms.min_samples:
            sample_dist = int(np.ceil(ref_len/self.prms.min_samples))
        elif approx_samps > self.prms.max_samples:
            sample_dist = int(np.floor(ref_len/self.prms.max_samples))
        else:
            sample_dist = self.prms.max_sample_dist

        fmlens = unc.self_align(self.prms.index_prefix, sample_dist)
        path_kfmlens = [p[self.prms.kmer_len-1:] if len(p) >= self.prms.kmer_len else [1] for p in fmlens]

        max_pathlen = 0
        all_pathlens = [len(p) for p in path_kfmlens if len(p) <= self.prms.max_replen]
        gt1_counts = np.zeros(max(all_pathlens))
        for l in all_pathlens:
            for i in range(l):
                gt1_counts[i] += 1

        max_pathlen = np.flatnonzero(gt1_counts / len(all_pathlens) <= self.prms.pathlen_percentile)[0]
        max_fmexp = int(np.log2(max([p[0] for p in path_kfmlens])))+1
        fm_path_mat = np.zeros((max_fmexp, max_pathlen))

        for p in path_kfmlens:
            for i in range(min(max_pathlen, len(p))):
                fm_path_mat[int(np.log2(p[i])), i] += 1
            for i in range(len(p), max_pathlen):
                fm_path_mat[0, i] += 1

        mean_fm_locs = list()
        for f in range(max_fmexp):
            loc_weights = fm_path_mat[f] / np.sum(fm_path_mat[f])
            mean_fm_locs.append(np.sum([i*loc_weights[i] for i in range(max_pathlen)]))
        self.fm_locs = np.array(mean_fm_locs)

        mean_loc_fms = list()
        for p in range(max_pathlen):
            fm_weights = fm_path_mat[:,p] / np.sum(fm_path_mat[:,p])
            mean_loc_fms.append(np.sum([i*fm_weights[i] for i in range(max_fmexp)]))
        self.loc_fms = np.array(mean_loc_fms)

        self.speed_denom = np.sum(self.loc_fms)

        self.prms_locs = np.arange(np.round(self.fm_locs[0]))
        self.all_locs = np.arange(max_pathlen)

    def get_model_threshs(self, fname=MODEL_THRESHS_FNAME):
        prob_thresh_in = open(fname)
        threshs = list()
        freqs = list()
        counts = list()
        for line in prob_thresh_in:
            thresh, freq, count = line.split()
            threshs.append(float(thresh))
            freqs.append(float(freq))
            counts.append(float(count))

        self.model_ekms = np.flip(np.array(threshs),0)
        self.model_pcks = np.flip(np.array(freqs),0)
        self.model_counts = np.flip(np.array(counts),0)

    def get_fn_speed(self, fn_locs, fn_pcks):
        pcks = np.interp(self.all_locs, fn_locs, fn_pcks)
        counts = np.interp(pcks, self.model_pcks, self.model_counts)
        speed = np.dot(counts, self.loc_fms) / (self.speed_denom)
        return speed

    def get_fn_prob(self, fn_locs, fn_pcks):
        return np.prod(np.interp(self.prms_locs, fn_locs, fn_pcks))

    def add_preset(self, name, tgt_prob=None, tgt_speed=None, exp_st=2, init_fac=2, eps=0.00001):

        exp = exp_st
        exp_min, exp_max = (None, None)

        pdelta = None

        pck1 = self.pck1
        pck2 = self.pck2

        sys.stderr.write("Computing %s parameters\n" % name)

        while True:
            fn_locs,fn_pcks = power_fn(self.fm_locs[0], pck1, pck2, exp)

            if tgt_prob is not None:
                delta = self.get_fn_prob(fn_locs, fn_pcks) - tgt_prob
            elif tgt_speed is not None:
                delta = self.get_fn_speed(fn_locs, fn_pcks) - tgt_speed

            if abs(delta) <= eps:
                break
            
            if delta == pdelta:
                #This works well for small references
                #TODO: check for larger references
                sys.stderr.write("Maxed out %s parameters\n" % name)
                break
            pdelta = delta

            if delta < 0:
                exp_max = exp
            else:
                exp_min = exp
            
            pexp = exp

            if exp_max == None:
                exp *= init_fac
            elif exp_min == None:
                exp /= init_fac
            else:
                exp = exp_min + ((exp_max - exp_min) / 2.0)

            #for floating point rounding errors
            if exp == pexp:
                break

        fm_pcks = np.interp(self.fm_locs, fn_locs, fn_pcks)
        fm_ekms = np.interp(fm_pcks, self.model_pcks, self.model_ekms)
        prob = self.get_fn_prob(fn_locs, fn_pcks)
        speed = self.get_fn_speed(fn_locs, fn_pcks)

        #while len(fm_ekms) > 2 and fm_ekms[-1] == fm_ekms[-2]:
        #    fm_ekms = fm_ekms[:-1]

        sys.stderr.write("Writing %s parameters\n" % name)
        self.functions[name] = (fm_ekms, prob, speed)

    def write(self):
        params_out = open(self.out_fname, "w")
        
        for name, fn in self.functions.items():
            ekms, prob, speed = fn
            params_out.write("%s\t%s\t%.5f\t%.3f\n" % (name, ",".join(map(str,ekms)), prob, speed))

        params_out.close()

