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

from _uncalled import _RefIndex, _RefCoord
from .argparse import Opt
from .config import ParamGroup
from . import nt

def str_to_coord(coord_str):
    spl = coord_str.split(":")         
    #name = spl[0]
    i = 1
    while True:
        try:
            st,en = map(int, spl[i].split("-"))
            break
        except ValueError:
            i += 1
    name = ":".join(spl[:i])
    
    if i < len(spl):
        return RefCoord(name,st,en)
    else:                              
        fwd = spl[i] == "+"
        return RefCoord(name,st,en,fwd)

class RefCoord(_RefCoord):
    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
            if isinstance(args, str):
                args = self._str_to_tuple(args)
            if not isinstance(args, tuple):
                raise ValueError(f"Invalid RefCoord: {args}")
        _RefCoord.__init__(self, *args)

    def _str_to_tuple(self, coord_str):
        spl = coord_str.split(":")         
        name = spl[0]                        
        st,en = map(int, spl[1].split("-"))
        
        if len(spl) == 2:                  
            return (name,st,en)
        elif len(spl) == 3 and spl[2] in {"+","-"}:
            fwd = spl[2] == "+"
            return (name,st,en,fwd)
        return None

class CoordSpace:

    def __init__(self, ref_name, refs, mrefs, fwd, kmers):
        self.ref_name = ref_name
        self.refs = refs
        self.mrefs = mrefs
        self.fwd = fwd

        self.set_kmers(kmers)

        if self.stranded and not isinstance(mrefs, pd.Index):
            raise ValueError("mrefs must be pandas.Index for stranded CoordSpace")
        if not self.stranded and not isinstance(mrefs, tuple):
            raise ValueError("mrefs must be tuple for stranded CoordSpace")

    def set_kmers(self, kmers):
        self.kmers = kmers
        if kmers is None:
            self.ref_kmers = None
        elif self.stranded:
            self.ref_kmers = pd.concat(
                {int(self.fwd) : self.kmers.set_axis(self.mref_to_ref(self.kmers.index))})
        else:
            self.ref_kmers = pd.concat(
                {0 : self.kmers[0].set_axis(self.mref_to_ref(self.kmers[0].index)), 
                 1 : self.kmers[1].set_axis(self.mref_to_ref(self.kmers[1].index))}
            )

    @property
    def stranded(self):
        return self.fwd is not None

    def contains(self, other):
        if isinstance(other, CoordSpace):
            return (self.ref_name == other.ref_name and
                    self.refs.min() <= other.refs.min() and 
                    self.refs.max() >= other.refs.max())
        elif isinstance(other, RefCoord):
            return (self.ref_name == other.name and
                    self.refs.min() <= other.start and 
                    self.refs.max() > other.end)
        else:
            raise ValueError(f"Expected CoordSpace or RefCoord, received {type(other)}")

    def _minmax_intersect(self, a, b):
        insc = a.intersection(b)
        return insc.min(),insc.max()

    def ref_slice(self, ref_start, ref_end):
        st = self.refs.get_loc(ref_start)
        en = self.refs.get_loc(ref_end-1)+1
        return self._slice(st,en)

    def mref_intersect(self, mrefs):
        if self.stranded:
            fwd_in = self.fwd
            minmax = self._minmax_intersect(self.mrefs, mrefs)
            if np.any(np.isnan(minmax)):
                return None
            bounds = self.mrefs.get_indexer(minmax)
        else:
            fwd_in = None
            for i in range(2):
                minmax = self._minmax_intersect(self.mrefs[i], mrefs)
                if not np.any(np.isnan(minmax)):
                    fwd_in = i
                    break
            if fwd_in is None:
                return None
            bounds = self.mrefs[fwd_in].get_indexer(minmax)

        st,en = sorted(bounds)
        en += 1

        return self._slice(st,en,fwd_in)

    def _slice(self, st, en, fwd=None):
        refs = self.refs[st:en]

        kmers = None

        if self.stranded:
            fwd = self.fwd
            mrefs = self.mrefs[st:en]
            if self.kmers is not None:
                kmers = self.kmers[st:en]

        elif fwd is None:
            mrefs = tuple( (m[st:en] for m in self.mrefs) )
            kmers = tuple( (k[st:en] for k in self.kmers) )
        else:
            mrefs = self.mrefs[fwd][st:en]
            if self.kmers is not None:
                kmers = self.kmers[fwd][st:en]

        #elif not stranded_out:

        return CoordSpace(self.ref_name, refs, mrefs, fwd, kmers)

    def ref_to_mref(self, ref, fwd=None):
        if isinstance(ref, (collections.abc.Sequence, np.ndarray, pd.Index)):
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

    #TODO check fwd and rev, return None if neither
    def is_mref_fwd(self, mref):
        if isinstance(mref, (collections.abc.Sequence, np.ndarray, pd.Index)):
            mref = mref[0]

        if self.stranded:
            return self.fwd and mref in self.mrefs
        else:
            return mref in self.mrefs[1]
            
    #TODO make single private method with fwd param. Also probably merge with above
            
    def all_mrefs_fwd(self, mrefs):

        if isinstance(mrefs, (collections.abc.Sequence, np.ndarray, pd.Index, pd.Series)):
            isin = lambda a,b: a.isin(b).all()
        else:
            isin = lambda a,b: a in b

        if self.stranded:
            if isin(mrefs, self.mrefs):
                return self.fwd
            return None
        for fwd,s_mrefs in enumerate(self.mrefs):
            if isin(mrefs, s_mrefs):
                return bool(fwd)
        return None
    
    def mref_to_ref(self, mref):
        fwd = self.all_mrefs_fwd(mref)

        if fwd is None:
            raise ValueError("mref coordinates outside of CoordSpace")

        if self.stranded:
            mrefs = self.mrefs
        else:
            mrefs = self.mrefs[fwd]

        if isinstance(mref, (collections.abc.Sequence, np.ndarray, pd.Index, pd.Series)):
            i = mrefs.get_indexer(mref)
        else:
            i = mrefs.get_loc(mref)

        return self.refs[i]

    def mref_to_ref_index(self, mrefs, multi=False):
        if self.stranded:
            refs = self.mref_to_ref(mrefs)
            fwd_mask = np.full(len(refs), self.fwd)

        else:
            revs = self.mrefs[0].get_indexer(mrefs)
            fwds = self.mrefs[1].get_indexer(mrefs)
            rev_mask = revs >= 0
            fwd_mask = fwds >= 0
            if not np.all(fwd_mask | rev_mask):
                raise ValueError("mrefs coordinates are outside of the CoordSpace")

            refs = np.zeros(len(mrefs), dtype=int)
            refs[rev_mask] = self.refs[revs[rev_mask]]
            refs[fwd_mask] = self.refs[fwds[fwd_mask]]

        if multi:
            strands = np.full(len(refs), "-")
            strands[fwd_mask] = "+"
            df = pd.DataFrame({
                "ref_name" : self.ref_name, 
                "ref" : refs,
                "strand" : strands})
            return pd.MultiIndex.from_frame(df)
        return pd.Index(refs, name="ref")

        #strand,refs = self.mref_to_ref(mrefs)
        #if str_strand:
        #    strand_label = "strand"
        #    strand = "+" if strand else "-"
        #else:
        #    strand_label = "fwd"
        #ret = pd.MultiIndex.from_product(
        #    [[self.ref_name], [, [strand]], 
        #    names=["ref_name", "ref", strand_label])

    #def validate_refs(self, refs):
    #    return len(self.refs.intersection(refs)) == len(refs)

    #def validate_mrefs(self, mrefs):
    #    return all_mrefs_fwd(mrefs) or all_mrefs_rev(mrefs)

    def __len__(self):
        return len(self.refs)

    def __repr__(self):
        return ("%s:%d-%d " % (self.ref_name, self.refs.start, self.refs.stop)) + str(self.mrefs)

class RefIndex(_RefIndex):

    def mrefs_to_kmers(self, mrefs, is_rna, kmer_trim):
        #if (mrefs.step < 0) == is_rna:
        #    #ref_coord.end -= kmer_shift
        #    kmers = self.get_kmers(mrefs.min(), mrefs.max()+1, is_rna)
        #else:
        #   #ref_coord.start += kmer_shift


        if kmer_trim:
            kmers = self.get_kmers(mrefs.min(), mrefs.max()+1, is_rna)
            if mrefs.step < 0:
                kmers = kmers[::-1]
            ret = pd.Series(index=mrefs[2:-2], data=kmers, name="kmer")
        else:
            kmers = self.get_kmers(mrefs.min()-2, mrefs.max()+3, is_rna)
            if mrefs.step < 0:
                kmers = kmers[::-1]
            ret = pd.Series(index=mrefs, data=kmers, name="kmer")
        return ret

    def ref_coord_to_mrefs(self, ref_coord, is_rna, flip_rev=True):
        st,en = self.ref_to_mref(ref_coord.name, ref_coord.start, ref_coord.end, ref_coord.fwd, is_rna)
        mrefs = pd.RangeIndex(st,en)
        if flip_rev and self.is_mref_flipped(mrefs.start):
            return mrefs[::-1]
        return mrefs

    #def full_coord_space(self, rid):

    #TODO get of -nt.K+1?
    def get_coord_space(self, ref_coord, is_rna, load_kmers=True, kmer_trim=False):

        rid = self.get_ref_id(ref_coord.name)
        if rid < 0:
            raise ValueError("Sequence not in reference: " + ref_coord.name)

        #if kmer_trim:
        #    ref_coord = RefCoord(ref_coord)
        #    ref_coord.start += 2
        #    ref_coord.end -= 2
        #if ref_coord.fwd is None or ref_coord.fwd == is_rna:
        #if is_rna:
        #    ref_coord.end -= kmer_shift
        #elif kmer_shift > 0:
        #    ref_coord.start += kmer_shift+1

        length = self.get_ref_len(rid)
        if ref_coord.start < 0 or ref_coord.end > length:
            raise ValueError("Reference coordinates %s out of bounds for sequence of length %d" % (ref_coord, length))

        refs = pd.RangeIndex(ref_coord.start, ref_coord.end)

        if ref_coord.stranded:
            mrefs = self.ref_coord_to_mrefs(ref_coord, is_rna)
        else:
            mrefs = tuple((
                self.ref_coord_to_mrefs(RefCoord(ref_coord, fwd), is_rna)
                for fwd in range(2)
            ))

        if load_kmers:
            if ref_coord.stranded:
                kmers = self.mrefs_to_kmers(mrefs, is_rna, kmer_trim)
            else:
                kmers = tuple(( 
                    self.mrefs_to_kmers(m, is_rna, kmer_trim) for m in mrefs
                ))
        else:
            kmers = None


        fwd = ref_coord.fwd if ref_coord.stranded else None 


        return CoordSpace(ref_coord.name, refs, mrefs, fwd, kmers)

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

