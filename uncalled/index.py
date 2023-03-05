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

import pandas as pd
from . import RefCoord, str_to_coord

import _uncalled
from .pore_model import PoreModel
from .argparse import Opt


class CoordSpace:

    def __init__(self, ref_name, refs, mposs, fwd, kmers, pac_shift):
        self.ref_name = ref_name
        self.refs = refs
        self.mposs = mposs
        self.fwd = fwd
        self.pac_shift = pac_shift
        self.pacs = refs + pac_shift

        self.set_kmers(kmers)

        if self.stranded:
            #if mposs.step < 0:
            #    self.refs = self.refs[::-1]
            #    self.pacs = self.pacs[::-1]
            if not isinstance(mposs, pd.Index):
                raise ValueError("mposs must be pandas.Index for stranded CoordSpace")
        if not self.stranded and not isinstance(mposs, tuple):
            raise ValueError("mposs must be tuple for stranded CoordSpace")

    def pac_to_pos(self, pac):
        return pac - self.pac_shift

    def pos_to_pac(self, ref):
        return ref + self.pac_shift

    def mpos_to_pac(self, mpos):
        return self.pac_shift + self.mpos_to_pos(mpos)

    def set_kmers(self, kmers):
        self.kmers = kmers
        if kmers is None:
            self.ref_kmers = None
        elif self.stranded:
            self.ref_kmers = pd.concat(
                {int(self.fwd) : self.kmers.set_axis(self.mpos_to_pos(self.kmers.index, self.fwd))})
        else:
            self.ref_kmers = pd.concat(
                {0 : self.kmers[0].set_axis(self.mpos_to_pos(self.kmers[0].index, False)), 
                 1 : self.kmers[1].set_axis(self.mpos_to_pos(self.kmers[1].index, True))}
            )

    @property
    def stranded(self):
        return self.fwd is not None

    @property
    def strands(self):
        return [int(self.fwd)] if self.stranded else [0,1]

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

    def ref_slice(self, ref_start=None, ref_end=None, fwd=None):
        if ref_start is not None:
            st = self.refs.get_loc(ref_start)
        else:
            st = 0
        if ref_end is not None:
            en = self.refs.get_loc(ref_end-1)+1
        else:
            en = len(self.refs)
        return self._slice(st,en,fwd)

    def mpos_intersect(self, mposs):
        if self.stranded:
            fwd_in = self.fwd
            minmax = self._minmax_intersect(self.mposs, mposs)
            if np.any(np.isnan(minmax)):
                return None
            bounds = self.mposs.get_indexer(minmax)
        else:
            fwd_in = None
            for i in range(2):
                minmax = self._minmax_intersect(self.mposs[i], mposs)
                if not np.any(np.isnan(minmax)):
                    fwd_in = i
                    break
            if fwd_in is None:
                return None
            bounds = self.mposs[fwd_in].get_indexer(minmax)

        st,en = sorted(bounds)
        en += 1

        return self._slice(st,en,fwd_in)

    def pac_intersect(self, coords):
        if self.stranded:
            try:
                pacs = coords[coords.get_loc(int(self.fwd))].get_level_values(1)
            except:
                return None
        else:
            pacs = coords.get_level_values(1)

        minmax = self._minmax_intersect(self.pacs, pacs)
        if np.any(np.isnan(minmax)):
            return None
        bounds = self.pacs.get_indexer(minmax)

        fwds = np.array(coords.get_level_values(0))
        all_fwd = np.all(fwds)
        all_rev = np.all(~fwds)
        fwd = None if not (all_fwd or all_rev) else all_fwd

        st,en = sorted(bounds)
        en += 1

        return self._slice(st,en,fwd)

    def _slice(self, st, en, fwd=None):
        refs = self.refs[st:en]

        kmers = None

        if self.stranded:
            fwd = self.fwd
            mposs = self.mposs[st:en]
            if self.kmers is not None:
                kmers = self.kmers[st:en]

        elif fwd is None:
            mposs = tuple( (m[st:en] for m in self.mposs) )
            if self.kmers is not None:
                kmers = tuple( (k[st:en] for k in self.kmers) )
        else:
            mposs = self.mposs[fwd][st:en]
            if self.kmers is not None:
                kmers = self.kmers[fwd].loc[st:en]

        #elif not stranded_out:

        return CoordSpace(self.ref_name, refs, mposs, fwd, kmers, self.pac_shift)

    def pos_to_mpos(self, ref, fwd=None):
        if isinstance(ref, (collections.abc.Sequence, np.ndarray, pd.Index)):
            i = self.refs.get_indexer(ref)
        else:
            i = self.refs.get_loc(ref)

        if self.stranded:
            if fwd is not None and fwd != self.fwd:
                raise ValueError("Invalid 'fwd' value or stranded CoordSpace")
            return self.mposs[i]

        # not stranded
        if fwd is None:
            return tuple( (mposs[i] for mposs in self.mposs) )
        return self.mposs[fwd][i]

    #TODO check fwd and rev, return None if neither
    def is_mpos_fwd(self, mpos):
        if isinstance(mpos, (collections.abc.Sequence, np.ndarray, pd.Index)):
            mpos = mpos[0]

        if self.stranded:
            return self.fwd and mpos in self.mposs
        else:
            return mpos in self.mposs[1]
            
    #TODO make single private method with fwd param. Also probably merge with above
            
    def all_mposs_fwd(self, mposs):

        if isinstance(mposs, (collections.abc.Sequence, np.ndarray, pd.Index, pd.Series)):
            isin = lambda a,b: a.dropna().isin(b).all()
        else:
            isin = lambda a,b: a in b

        if self.stranded:
            if isin(mposs, self.mposs):
                return self.fwd
            return None
        for fwd,s_mposs in enumerate(self.mposs):
            if isin(mposs, s_mposs):
                return bool(fwd)
        return None
    
    def mpos_to_pos(self, mpos, fwd=None):
        fwd = fwd if fwd is not None else self.all_mposs_fwd(mpos)

        if fwd is None:
            raise ValueError("mpos coordinates outside of CoordSpace")

        if self.stranded:
            mposs = self.mposs
        else:
            mposs = self.mposs[fwd]

        if isinstance(mpos, (collections.abc.Sequence, np.ndarray, pd.Index, pd.Series)):
            i = mposs.get_indexer(mpos)
        else:
            i = mposs.get_loc(mpos)

        return self.refs[i]

    def mpos_to_pos_index(self, mposs, multi=False):
        if self.stranded:
            refs = self.mpos_to_pos(mposs, self.fwd)
            fwd_mask = np.full(len(refs), self.fwd)

        else:
            revs = self.mposs[0].get_indexer(mposs)
            fwds = self.mposs[1].get_indexer(mposs)
            rev_mask = revs >= 0
            fwd_mask = fwds >= 0
            if not np.all(fwd_mask | rev_mask):
                raise ValueError("mposs coordinates are outside of the CoordSpace")

            refs = np.zeros(len(mposs), dtype=int)
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

        #strand,refs = self.mpos_to_pos(mposs)
        #if str_strand:
        #    strand_label = "strand"
        #    strand = "+" if strand else "-"
        #else:
        #    strand_label = "seq.fwd"
        #ret = pd.MultiIndex.from_product(
        #    [[self.ref_name], [, [strand]], 
        #    names=["ref_name", "ref", strand_label])

    #def validate_refs(self, refs):
    #    return len(self.refs.intersection(refs)) == len(refs)

    #def validate_mposs(self, mposs):
    #    return all_mposs_fwd(mposs) or all_mposs_rev(mposs)

    def __len__(self):
        return len(self.refs)

    def __repr__(self):
        return ("%s:%d-%d " % (self.ref_name, self.refs.start, self.refs.stop)) + str(self.mposs)

class RefIndex:

    def __init__(self, k, *args, **kwargs):

        self.InstanceClass = getattr(_uncalled, f"RefIndexK{k}", None)
        self.Model = getattr(_uncalled, f"PoreModelK{k}", None)
        if self.InstanceClass is None or self.Model is None:
            raise ValueError(f"Invalid k-mer length {k}")
        
        if "kmer_shift" in kwargs:
            shift = kwargs["kmer_shift"]
            del kwargs["kmer_shift"]
        else:
            shift = PoreModel.get_kmer_shift(k)

        self.trim = (shift, k-shift-1)

        self.instance = self.InstanceClass(*args, **kwargs)

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def mposs_to_kmers(self, ref_id, mposs, is_rna, kmer_trim):
        #if (mposs.step < 0) == is_rna:
        #    #ref_coord.end -= kmer_shift
        #    kmers = self.get_kmers(mposs.min(), mposs.max()+1, is_rna)
        #else:
        #   #ref_coord.start += kmer_shift

        if (mposs.step > 0) != is_rna:
            i,j = self.trim
        else:
            i,j = reversed(self.trim)

        if kmer_trim:
            kmers = self.get_kmers(ref_id, mposs.min(), mposs.max()+1, is_rna)
            if mposs.step < 0:
                kmers = kmers[::-1]
            ret = pd.Series(index=mposs[i:-j], data=kmers, name="kmer")
        else:
            kmers = self.get_kmers(ref_id, mposs.min()-i, mposs.max()+j+1, is_rna)
            if mposs.step < 0:
                kmers = kmers[::-1]
            ret = pd.Series(index=mposs, data=kmers, name="kmer")
        return ret

    def ref_coord_to_mposs(self, ref_coord, is_rna, flip_rev=True):
        st,en = self.pos_to_mpos(ref_coord.name, ref_coord.start, ref_coord.end, ref_coord.fwd, is_rna)
        mposs = pd.RangeIndex(st,en)
        if flip_rev and self.is_mpos_flipped(mposs.start):
            return mposs[::-1]
        return mposs

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
            mposs = self.ref_coord_to_mposs(ref_coord, is_rna)
        else:
            mposs = tuple((
                self.ref_coord_to_mposs(RefCoord(ref_coord, fwd), is_rna)
                for fwd in [False, True]
            ))

        if load_kmers:
            if ref_coord.stranded:
                kmers = self.mposs_to_kmers(rid, mposs, is_rna, kmer_trim)
            else:
                kmers = tuple(( 
                    self.mposs_to_kmers(rid, m, is_rna, kmer_trim) for m in mposs
                ))
        else:
            kmers = None


        fwd = ref_coord.fwd if ref_coord.stranded else None 

        pac_shift = self.get_pac_offset(ref_coord.name)

        return CoordSpace(ref_coord.name, refs, mposs, fwd, kmers, pac_shift)

_index_cache = dict()

def load_index(k, prefix, load_pacseq=True, load_bwt=False, cache=True, kmer_shift=None):
    idx = _index_cache.get(prefix, None)
    if idx is None:
        idx = RefIndex(k, prefix, load_pacseq, load_bwt, kmer_shift=kmer_shift)
        if cache: _index_cache[prefix] = idx
    else:
        if load_pacseq and not idx.pacseq_loaded():
            idx.load_pacseq()
        if load_bwt and not idx.bwt_loaded():
            idx.load_index(prefix)
    return idx



def index(conf):
    """Build an index from a FASTA reference"""
        
#       Uses BWA to generate pacseq and FM index files, then computes reference-specific probability thresholds. The FM index and probability thresholds are only required for Real-Time Enrichment subcommands (see --no-bwt if only interested in DTW analysis)"""

    prms = conf.index

    if prms.index_prefix is None or len(prms.index_prefix) == 0:
        prms.index_prefix = prms.fasta_filename

    bwa_built = True

    for suff in UNCL_SUFFS:
        if not os.path.exists(prms.index_prefix + suff):
            bwa_built = False
            break

    if bwa_built:
        sys.stderr.write("Using previously built BWA index.\nNote: to fully re-build the index delete files with the \"%s.*\" prefix.\n" % prms.index_prefix)
    else:
        _uncalled.RefIndexK5.create(prms.fasta_filename, prms.index_prefix, prms.no_bwt)
        
        if prms.no_bwt: 
            sys.stderr.write("Pacseq built\n")
            return

    sys.stderr.write("Initializing parameter search\n")
    p = IndexParameterizer(prms)

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

        fmlens = _uncalled.self_align(self.prms.index_prefix, sample_dist)
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

