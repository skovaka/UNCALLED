"""Compute, compare, and query alignment layers

subcommand options:
compare  Compare different alignment methods on the same set of reads"""

import sys, os
import numpy as np
import argparse
import re
import time
import scipy.stats.mstats as mstats
import types
import pandas as pd
import scipy.stats

from .. import config
from ..fast5 import parse_read_ids
from ..dtw import Tracks
from ..dtw.layers import parse_layers


def compare(conf):
    """Compute distance between alignments of the same reads"""
    t = time.time()

    group_b = "bcaln" if conf.bcaln else "dtw"

    if conf.bcaln:
        conf.tracks.layers += [("bc_cmp", "mean_ref_dist")] 
    else:
        conf.tracks.layers += [("cmp", "mean_ref_dist")] 

    all_layers = not (conf.jaccard or conf.mean_ref_dist)
    calc_jaccard = all_layers or conf.jaccard
    calc_mean_ref_dist = all_layers or conf.mean_ref_dist

    tracks = Tracks(conf=conf)

    t = time.time()
    t_all = time.time()

    t = time.time()

    for read_id,chunk in tracks.iter_reads():
        if chunk.any_empty:
            sys.stderr.write(f"Skipping {read_id}\n")
        else:
            #if conf.save:
            #    print(read_id)
            chunk.calc_compare(group_b, calc_jaccard, calc_mean_ref_dist, conf.save)
            chunk.write_alignment()

        print(f"{read_id}\t{time.time()-t:.4f}")
        sys.stdout.flush()
        t = time.time()

def dump(conf):
    """Output DTW alignment paths and statistics"""

    tracks = Tracks(conf=conf)
    #TODO add layer dependencies (compare requires start/length)
    #tracks.set_layers(["start", "length", "bcaln.start", "bcaln.length"] + conf.layers)
    tracks.set_layers(conf.layers)

    layer_groups = {group for group,_ in parse_layers(conf.layers)}

    need_cmp = "cmp" in layer_groups
    need_bc_cmp = "bc_cmp" in layer_groups

    header = True


    for read_id,tracks in tracks.iter_reads():
        for track in tracks:
            if header:
                columns = track.layers.index.names + [
                    ".".join([c for c in col if len(c) > 0]) 
                    for col in track.layers.columns]
                print("\t".join(columns))
                header = False
            layers = track.layers.dropna(axis=0, how="any")
            sys.stdout.write(layers.to_csv(sep="\t", header=False))
        t = time.time()
