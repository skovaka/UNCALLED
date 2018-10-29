#!/usr/bin/env python
import sys                         
import os
import numpy as np
import argparse
from bisect import bisect_left, bisect_right

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
MODEL_FNAME = os.path.join(ROOT_DIR, "models/r94_5mers.txt")
MODEL_THRESHS_FNAME = os.path.join(ROOT_DIR, "models/r94_5mers_threshs.txt")
PARAM_SUFF = ".uncl"


def get_fmlen_percentiles(fmlens, perc, max_path_len):
    path_fmlens = [list() for i in range(max_path_len)]
    path_percs = list()

    for path_lens in fmlens:
        for i in range(min(max_path_len, len(path_lens))):
            path_fmlens[i].append(path_lens[i])

        for i in range(len(path_lens), max_path_len):
            path_fmlens[i].append(1)

    for p in range(len(path_fmlens)):
        path_fmlens[p] = list(sorted(path_fmlens[p]))
        i = int(len(path_fmlens[p])*(perc/100))
        path_percs.append(path_fmlens[p][i])

        if path_percs[-1] == 1:
            break

    return path_percs

def power_fn(N, st, en, ex):
    return ((np.arange(N)/(N-1))**ex) * (en-st) + st

def match_prfn_em(fn, N, st, en, tgt_prod, eps=0.00001, param_st=2, init_fac=2):
    
    param = param_st
    param_min, param_max = (None, None)

    delta = np.prod(fn(N, st, en, param)) - tgt_prod

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

        delta = np.prod(fn(N, st, en, param)) - tgt_prod

    return param

def get_params(bwa_prefix, fmlens, kmer_len, 
                 fm_percentile, matchpr1, matchpr2, match_prod):


    path_fmlens = get_fmlen_percentiles(fmlens,
                                   fm_percentile, 25)[kmer_len-1:]
    path_len = len(path_fmlens)

    prob_thresh_in = open(MODEL_THRESHS_FNAME)
    prob_freqs = dict()
    freq_probs = list()
    for line in prob_thresh_in:
        prob, match_freq, kmer_count = line.split()
        prob_freqs[float(prob)] = float(match_freq)
        freq_probs.append((float(match_freq), float(prob)))

    freq_probs = list(sorted(freq_probs))

    exp = match_prfn_em(power_fn, 
                        path_len, 
                        matchpr1, 
                        matchpr2, 
                        match_prod)

    fn_freqs = power_fn(path_len, matchpr1, matchpr2, exp)

    out_fname = bwa_prefix + PARAM_SUFF
    params_out = open(out_fname, "w")
    sys.stderr.write("Writing %s\n" % out_fname)
    
    for i in range(path_len):
        if path_fmlens[i] > 1 and path_fmlens[i] == path_fmlens[i+1]:
            continue
        j = bisect_left(freq_probs, (fn_freqs[i], 0))
        if i > 0:
            params_out.write("%d\t" % (path_fmlens[i]))
        params_out.write("%.2f\t%.4f\n" % (freq_probs[j][1], fn_freqs[i]))

