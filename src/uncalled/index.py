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
import argparse
from uncalled import mapping, params
from bisect import bisect_left, bisect_right


def power_fn(xmax, ymin, ymax, exp, N=100):
    dt = 1.0/N
    t = np.arange(0, 1+dt, dt)

    return t*xmax, (t**exp) * (ymax-ymin) + ymin

class IndexParameterizer:

    def __init__(self, args):
        self.out_fname = args.bwa_prefix + params.INDEX_UNCL_SUFF

        self.pck1 = args.matchpr1
        self.pck2 = args.matchpr2

        self.calc_map_stats(args)
        self.get_model_threshs()

        self.functions = dict()

    def calc_map_stats(self, args):

        ann_in = open(args.bwa_prefix + params.INDEX_ANN_SUFF)
        header = ann_in.readline()
        ref_len = int(header.split()[0])
        ann_in.close()

        approx_samps = ref_len / args.max_sample_dist
        if approx_samps < args.min_samples:
            sample_dist = int(np.ceil(ref_len/args.min_samples))
            sys.stderr.write("Maxed %.2f\n" % sample_dist)
        elif approx_samps > args.max_samples:
            sample_dist = int(np.floor(ref_len/args.max_samples))
        else:
            sample_dist = args.max_sample_dist
            sys.stderr.write("NOT maxed %.2f\n" % sample_dist)

        fmlens = mapping.self_align(args.bwa_prefix, sample_dist)
        path_kfmlens = [p[args.kmer_len-1:] if len(p) >= args.kmer_len else [1] for p in fmlens]

        max_pathlen = 0
        all_pathlens = [len(p) for p in path_kfmlens if len(p) <= args.max_replen]
        gt1_counts = np.zeros(max(all_pathlens))
        for l in all_pathlens:
            for i in range(l):
                gt1_counts[i] += 1

        max_pathlen = np.flatnonzero(gt1_counts / len(all_pathlens) <= args.pathlen_percentile)[0]
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

        self.conf_locs = np.arange(np.round(self.fm_locs[0]))
        self.all_locs = np.arange(max_pathlen)

    def get_model_threshs(self, fname=params.MODEL_THRESHS_FNAME):
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
        return np.prod(np.interp(self.conf_locs, fn_locs, fn_pcks))

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

