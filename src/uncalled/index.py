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
from bisect import bisect_left, bisect_right

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
MODEL_FNAME = os.path.join(ROOT_DIR, "models/r94_5mers.txt")
MODEL_THRESHS_FNAME = os.path.join(ROOT_DIR, "models/r94_5mers_threshs.txt")
PARAM_SUFF = ".uncl"

def get_model_threshs(fname=MODEL_THRESHS_FNAME):
    prob_thresh_in = open(fname)
    threshs = list()
    freqs = list()
    counts = list()
    for line in prob_thresh_in:
        thresh, freq, count = line.split()
        threshs.append(float(thresh))
        freqs.append(float(freq))
        counts.append(float(count))

    return (list(reversed(threshs)),
            list(reversed(freqs)),
            list(reversed(counts)))

def get_max_pathlen(path_fmlens, pathlen_perc):
    pathlens = [len(p) for p in path_fmlens]

    gt1_counts = np.zeros(max(pathlens))
    for l in pathlens:
        for i in range(l):
            gt1_counts[i] += 1
    gt1_fracs = gt1_counts / len(pathlens)

    max_pathlen = 0
    while gt1_counts[max_pathlen] / len(pathlens) > pathlen_perc:
        max_pathlen += 1
        if max_pathlen >= len(gt1_counts):
            sys.stderr.write("Error: reference too repetitive for option -m=%f\n" % args.max_multi_frac)
            sys.exit(1)

    return max_pathlen

def get_mean_fm_locs(path_kfmlens, pathlen_perc):

    max_pathlen = get_max_pathlen(path_kfmlens, pathlen_perc)

    max_fmexp = int(np.log2(max([p[0] for p in path_kfmlens])))+1


    fm_path_mat = np.zeros((max_fmexp, max_pathlen))

    for p in path_kfmlens:
        for i in range(min(max_pathlen, len(p))):
            fm_path_mat[int(np.log2(p[i])), i] += 1
        for i in range(len(p), max_pathlen):
            fm_path_mat[0, i] += 1

    mean_pathlocs = list()
    for f in range(max_fmexp):
        pathlen_freqs = fm_path_mat[f] / np.sum(fm_path_mat[f])
        mean_pathlocs.append(np.sum([(i)*pathlen_freqs[i] for i in range(max_pathlen)]))
    mean_pathlocs = np.array(mean_pathlocs)

    return mean_pathlocs

def power_fn(pathlen, pr1, pr2, exp, N=100):
    dt = 1.0/N
    t = np.arange(0, 1+dt, dt)

    return t*pathlen, (t**exp) * (pr2-pr1) + pr1

def matchfn_em(fn, pathlen, pr1, pr2, tgt_prod, 
                  param_st=2, init_fac=2, hyperparams=(), eps=0.00001):
   
    param = param_st
    param_min, param_max = (None, None)

    intpath = np.arange(np.round(pathlen))

    px,py = fn(pathlen, pr1, pr2, param, *hyperparams)
    delta = np.prod([np.interp(p, px, py) for p in intpath]) - tgt_prod

    while abs(delta) > eps:
        if delta < 0:
            param_max = param
        else:
            param_min = param

        if param_max == None:
            param *= init_fac
        elif param_min == None:
            param /= init_fac
        else:
            param = param_min + ((param_max - param_min) / 2.0)

        px,py = fn(pathlen, pr1, pr2, param, *hyperparams)
        #delta = np.prod([interpolate(p, px, py) for p in intpath]) - tgt_prod
        delta = np.prod([np.interp(p, px, py) for p in intpath]) - tgt_prod

    return param,px,py

def get_params(bwa_prefix, fmlens, kmer_len, 
                 fm_percentile, matchpr1, matchpr2, match_prod):

    path_kfmlens = [p[kmer_len-1:] for p in fmlens]


    model_ekms,model_pcks,model_counts = get_model_threshs()

    mean_fm_locs = get_mean_fm_locs(path_kfmlens, fm_percentile)

    pwr_exp,pwr_locs,pwr_pcks = matchfn_em(
                    power_fn, 
                    mean_fm_locs[0], 
                    matchpr1, 
                    matchpr2, 
                    match_prod)

    fm_pcks = np.interp(mean_fm_locs, pwr_locs, pwr_pcks)

    out_fname = bwa_prefix + PARAM_SUFF
    params_out = open(out_fname, "w")
    params_out.write("power\t%s\n" % ",".join(map(str,np.interp(fm_pcks, model_pcks, model_ekms))))
    params_out.close()
    #sys.stderr.write("Writing %s\n" % out_fname)

    #path_fmlens[-1]=1
    #
    #for i in range(path_len):
    #    if path_fmlens[i] > 1 and path_fmlens[i] == path_fmlens[i+1]:
    #        continue
    #    j = bisect_left(freq_probs, (fn_freqs[i], 0))
    #    if i > 0:
    #        params_out.write("%d\t" % (path_fmlens[i]))
    #    params_out.write("%.2f\t%.4f\n" % (freq_probs[j][1], fn_freqs[i]))

