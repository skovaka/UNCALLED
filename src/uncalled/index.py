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
from uncalled import mapping
from bisect import bisect_left, bisect_right
import matplotlib.pyplot as plt

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

    return (np.flip(np.array(threshs)),
            np.flip(np.array(freqs)),
            np.flip(np.array(counts)))

def get_conf_pathlen(path_fmlens, pathlen_perc, max_pathlen):
    pathlens = [len(p) for p in path_fmlens if len(p) <= max_pathlen]

    gt1_counts = np.zeros(max(pathlens))
    for l in pathlens:
        for i in range(l):
            gt1_counts[i] += 1
    gt1_fracs = gt1_counts / len(pathlens)

    conf_pathlen = 0
    while gt1_counts[conf_pathlen] / len(pathlens) > pathlen_perc:
        conf_pathlen += 1

    return conf_pathlen
    #return 12

def get_mean_fm_locs(path_kfmlens, pathlen_perc, max_pathlen):


    conf_pathlen = get_conf_pathlen(path_kfmlens, pathlen_perc, max_pathlen)

    max_fmexp = int(np.log2(max([p[0] for p in path_kfmlens])))+1


    fm_path_mat = np.zeros((max_fmexp, conf_pathlen))

    for p in path_kfmlens:
        for i in range(min(conf_pathlen, len(p))):
            fm_path_mat[int(np.log2(p[i])), i] += 1
        for i in range(len(p), conf_pathlen):
            fm_path_mat[0, i] += 1

    mean_fm_locs = list()
    for f in range(max_fmexp):
        pathlen_freqs = fm_path_mat[f] / np.sum(fm_path_mat[f])
        mean_fm_locs.append(np.sum([i*pathlen_freqs[i] for i in range(conf_pathlen)]))
    mean_fm_locs = np.array(mean_fm_locs)

    mean_loc_fms = list()
    for p in range(conf_pathlen):
        fmexp_freqs = fm_path_mat[:,p] / np.sum(fm_path_mat[:,p])
        mean_loc_fms.append(np.sum([i*fmexp_freqs[i] for i in range(max_fmexp)]))
    mean_loc_fms = np.array(mean_loc_fms)
    return mean_fm_locs, mean_loc_fms

def power_fn(pathlen, pr1, pr2, exp, N=100):
    dt = 1.0/N
    t = np.arange(0, 1+dt, dt)

    return t*pathlen, (t**exp) * (pr2-pr1) + pr1

def bezier(pts, N=100, t=None):
    if t is None:
        dt = 1.0/N
        t = np.arange(0, 1+dt, dt)

    if len(pts.shape) == 2:
        return bezier(pts[:,0],t=t), bezier(pts[:,1],t=t)

    if len(pts) == 1:
        return pts[0:]

    return (1-t) * bezier(pts[:-1],t=t) + t * bezier(pts[1:],t=t)

def bezier_quad(pathlen, pr1, pr2, dist, skew):

    if dist == 1:
        dx1=dx2=sx = pathlen
        dy1=dy2=sy = pr1
    else:
        dx1 = dist*pathlen
        dy1 = pr1
        dx2 = pathlen
        dy2 = pr1+(1-dist)*(pr2-pr1)
        sx = dx1 + skew * (dx2 - dx1)
        sy = dy1 + skew * (dy2 - dy1)

    bx, by = bezier(np.array([[0, pr1], [sx, sy], [pathlen, pr2]]))

    return bx,by

def bezier_cubic_cw(pathlen, pr1, pr2, center, width, plot_col=None):

    cx = pr1 + (pr2 - pr1)*center

    x1 = cx + width * (pathlen-cx)
    x2 = cx - width * cx

    bx, by = bezier(np.array([[0, pr1], [x1, pr1], [x2, pr2], [pathlen, pr2]]))
    
    if plot_col != None:
        plt.plot(bx,by, color=plot_col)
        plt.scatter([x1, x2], [pr1, pr2], c=plot_col)

    return bx,by

def linear(pathlen, pr1, pr2, prm):
    xs = np.arange(pathlen)
    ys = np.zeros(len(xs))

    ys[:-1] = np.interp(xs[:-1], [0, pathlen-1], [pr1, pr1+(pr2-pr1)*(1-prm)])
    ys[-1] = pr2

    return xs,ys

def matchfn_em(fn, pathlen, pr1, pr2, tgt_prod=None, 
               tgt_speed=None, model_pcks=None, model_counts=None,
               param_st=2, init_fac=2, hyperparams=(), eps=0.00001):
   
    param = param_st
    param_min, param_max = (None, None)

    intpath = np.arange(np.round(pathlen))

    fn_locs,fn_pcks = fn(pathlen, pr1, pr2, param, *hyperparams)
    path_pcks = np.interp(intpath, fn_locs, fn_pcks)

    if tgt_prod is not None:
        delta = np.prod(path_pcks) - tgt_prod
    elif tgt_speed is not None:
        path_counts = p.interp(path_pcks, model_pcks, model_counts)
        #speed_fac = np.sum([c*f for c,f in zip(path_counts, path_fms)]) / (np.arange(len(intpath))+1).sum()


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
        delta = np.prod([np.interp(p, px, py) for p in intpath]) - tgt_prod

    return param,px,py


def bezier_cubic_search(pathlen, pr1, pr2, tgt_prod, n, eps):
    fns = list()
    
    p0 = [0, pr1]
    p3 = [pathlen, pr2]

    dx = pathlen / n
    dy = (pr2 - pr1) / n

    intpath = np.arange(np.round(pathlen))

    for x1 in np.arange(0, pathlen+dx, dx):
        for y1 in np.arange(pr1, pr2+dy, dy):
            for x2 in np.arange(0, pathlen+dx, dx):
                for y2 in np.arange(pr1, pr2+dy, dy):
                    bx, by = bezier(np.array([p0,[x1,y1],[x2,y2],p3]))
                    delta = abs(np.prod([np.interp(p, bx, by) for p in intpath]) - tgt_prod)
                    if delta < eps:
                        fns.append( (p0[0],p0[1],x1,y1,x2,y2,p3[0],p3[1]) )

    return fns

def write_params(args):

    out_fname = args.bwa_prefix + PARAM_SUFF
    params_out = open(out_fname, "w")

    model_ekms,model_pcks,model_counts = get_model_threshs()

    fmlens = mapping.self_align(args.bwa_prefix, args.ref_fasta, args.sample_dist)
    path_kfmlens = [p[args.kmer_len-1:] if len(p) >= args.kmer_len else [1] for p in fmlens]
    mean_fm_locs, mean_loc_fms = get_mean_fm_locs(path_kfmlens, args.pathlen_percentile, args.max_pathlen)

    sys.stdout.write("Max path length: %f\n" % mean_fm_locs[0])

    prm,locs,pcks = matchfn_em(
                        power_fn, 
                        mean_fm_locs[0], 
                        args.matchpr1, 
                        args.matchpr2, 
                        args.target_prod)
    fm_pcks = np.interp(mean_fm_locs, locs, pcks)
    fm_ekms = np.interp(fm_pcks, model_pcks, model_ekms)
    params_out.write("power\t%s\n" % (",".join(map(str,fm_ekms))))

    intpath = np.arange(np.round(mean_fm_locs[0]))
    fm_exps = np.arange(len(mean_fm_locs))
    path_fms = np.interp(intpath, fm_exps, mean_fm_locs)
    #print(mean_fm_locs)

    for tgt in np.arange(0.02, 0.201, 0.005):
        prm,locs,pcks = matchfn_em(
                            power_fn, 
                            mean_fm_locs[0], 
                            args.matchpr1, 
                            args.matchpr2, 
                            tgt)

        fm_pcks = np.interp(mean_fm_locs, locs, pcks)
        fm_ekms = np.interp(fm_pcks, model_pcks, model_ekms)

        path_pcks = np.interp(intpath, locs, pcks)
        path_counts = np.interp(path_pcks, model_pcks, model_counts)

        speed_fac = np.sum([c*f for c,f in zip(path_counts, path_fms)]) / (np.arange(len(intpath))+1).sum()

        
        #TODO: compute this better
        path_counts2 = np.full(len(mean_loc_fms), path_counts[-1])
        path_counts2[:len(path_counts)] = path_counts

        mlf=mean_loc_fms
        #speed_fac2 = np.sum([c*f for c,f in zip(path_counts2, mlf)]) / (np.sum(mlf+1))
        speed_fac2 = np.dot(path_counts2, mlf) / (np.sum(mlf)+len(mlf))

        #speed_fac2 = np.sum([c*f for c,f in zip(path_counts2, mean_loc_fms)]) / (np.sum(mean_loc_fms))
        #print(speed_fac2)

        params_out.write("power_%.3f\t%s\t%.2f\n" % (tgt, ",".join(map(str,fm_ekms)), speed_fac2))
    
    plt.plot(locs,pcks)
    #plt.plot(mean_fm_locs,fm_ekms)
    #prm,locs,pcks = logprfn_em(model_ekms, model_pcks,
    #                    power_fn,
    #                    mean_fm_locs[0],
    #                    -10.07,
    #                    -2.268,
    #                    args.target_prod)
    #fm_ekms = np.flip(np.interp(mean_fm_locs, locs, pcks))
    #params_out.write("logpr_power\t%s\n" % (",".join(map(str,fm_ekms))))

    #fns = bezier_cubic_search(mean_fm_locs[0], 
    #                    args.matchpr1, 
    #                    args.matchpr2, 
    #                    args.target_prod, 
    #                    15, 
    #                    0.001)

    ##sys.stderr.write("%d functions found\n" % len(fns))

    #for x0,y0,x1,y1,x2,y2,x3,y3 in fns:
    #    locs,pcks = bezier(np.array([[x0,y0],[x1,y1],[x2,y2],[x3,y3]]))
    #    fm_pcks = np.interp(mean_fm_locs, locs, pcks)
    #    fm_ekms = np.interp(fm_pcks, model_pcks, model_ekms)
    #    params_out.write("search_%f,%f_%f,%f_%f,%f_%f,%f\t%s\n" % (x0,y0,x1,y1,x2,y2,x3,y3, ",".join(map(str,fm_ekms))))

    #for s in np.arange(0, 1.05, 0.05):
    #    #sys.stderr.write("Creating bezier quad with skew=%.2f -> " % s)
    #    exp,locs,pcks = matchfn_em(
    #                        bezier_quad, 
    #                        mean_fm_locs[0], 
    #                        args.matchpr1, 
    #                        args.matchpr2, 
    #                        args.target_prod,
    #                        0.5,
    #                        hyperparams=(s,))
    #    fm_pcks = np.interp(mean_fm_locs, locs, pcks)
    #    fm_ekms = np.interp(fm_pcks, model_pcks, model_ekms)
    #    params_out.write("bezier_q%.2f\t%s\n" % (s, ",".join(map(str,fm_ekms))))

    #for w in np.arange(-0.9, 1.0, 0.1):
    #    #sys.stderr.write("Creating bezier cubic with width=%.1f\n" % w)
    #    exp,locs,pcks = matchfn_em(
    #                        bezier_cubic_cw, 
    #                        mean_fm_locs[0], 
    #                        args.matchpr1, 
    #                        args.matchpr2, 
    #                        args.target_prod,
    #                        0.5,
    #                        hyperparams=(w,))
    #    fm_pcks = np.interp(mean_fm_locs, locs, pcks)
    #    fm_ekms = np.interp(fm_pcks, model_pcks, model_ekms)
    #    params_out.write("bezier_cw%.1f\t%s\n" % (w, ",".join(map(str,fm_ekms))))

    
   # exp,locs,pcks = matchfn_em(
   #                     bezier_cubic_cw, 
   #                     mean_fm_locs[0], 
   #                     args.matchpr1, 
   #                     args.matchpr2, 
   #                     args.target_prod,
   #                     0.5,
   #                     hyperparams=(0.9,))
   # bezier_cubic_cw(mean_fm_locs[0], args.matchpr1, args.matchpr2, exp, 0.9, 'red')

   # exp,locs,pcks = matchfn_em(
   #                     bezier_cubic_cw, 
   #                     mean_fm_locs[0], 
   #                     args.matchpr1, 
   #                     args.matchpr2, 
   #                     args.target_prod,
   #                     0.5,
   #                     hyperparams=(-0.6,))
   # bezier_cubic_cw(mean_fm_locs[0], args.matchpr1, args.matchpr2, exp, -0.6, 'blue')
   # fm_pcks = np.interp(mean_fm_locs, locs, pcks)
   # fm_ekms = np.interp(fm_pcks, model_pcks, model_ekms)
    #plt.plot(mean_fm_locs,fm_ekms)

    #plt.legend(['bez_cw0.9', 'bez_cw-0.6'])

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

