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
from ..dtw import Tracks
from ..dtw.layers import parse_layers


def compare(conf):
    """Compute distance between alignments of the same reads"""
    t = time.time()

    group_b = "moves" if conf.moves else "dtw"

    if conf.moves:
        conf.tracks.layers += [("mvcmp", "dist")] 
    else:
        conf.tracks.layers += [("cmp", "dist")] 

    all_layers = not (conf.jaccard or conf.dist)
    calc_jaccard = all_layers or conf.jaccard
    calc_dist = all_layers or conf.dist

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
            chunk.calc_compare(group_b, False, conf.save)
            chunk.write_alignment()

        print(f"{read_id}\t{time.time()-t:.4f}")
        sys.stdout.flush()
        t = time.time()

def dump(conf):
    """Output DTW alignment paths and statistics"""

    tracks = Tracks(conf=conf)
    #TODO add layer dependencies (compare requires start/length)
    #tracks.set_layers(["start", "length", "moves.start", "moves.length"] + conf.layers)
    tracks.set_layers(conf.layers)

    layer_groups = {group for group,_ in parse_layers(conf.layers)}

    need_cmp = "cmp" in layer_groups
    need_mvcmp = "mvcmp" in layer_groups

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
