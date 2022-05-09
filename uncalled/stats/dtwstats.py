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
from ..sigproc import ProcRead
from ..argparse import Opt, comma_split, ref_coords
from ..index import BWA_OPTS
from ..fast5 import parse_read_ids
from ..dtw import Tracks, LAYERS
from ..dtw.aln_track import parse_layers

COMPARE_OPTS = (
    Opt("db_in", "tracks.io"),
    Opt(("-l", "--read-filter"), "tracks", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "tracks"),
    Opt(("-b", "--bcaln"), action="store_true", help="Compare against basecalled alignment. If two tracks input will look for \"bcaln\" group in second track, otherwise will look in the first track."),
    Opt(("-s", "--save"), action="store_true", help="Will save in database if included, otherwise will output to TSV"),
    Opt(("-j", "--jaccard"), action="store_true", help="Will compute per-reference raw sample jaccard distances. Output by default if no other statistics are specified."),
    Opt(("-d", "--mean-ref-dist"), action="store_true", help="Will compute mean reference coordinate distances between raw samples of alignments of the same read. Output by default if no other statistics are specified."),
    #Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
)

def compare(conf):
    """Compute distance between alignments of the same reads"""
    t = time.time()

    group_b = "bcaln" if conf.bcaln else "dtw"

    if conf.bcaln:
        conf.tracks.layers += [("bc_cmp", "mean_ref_dist")] #LAYERS["bc_cmp"]["mean_ref_dist"].deps
    else:
        conf.tracks.layers += [("cmp", "mean_ref_dist")] #LAYERS["cmp"]["mean_ref_dist"].deps

    all_layers = not (conf.jaccard or conf.mean_ref_dist)
    calc_jaccard = all_layers or conf.jaccard
    calc_mean_ref_dist = all_layers or conf.mean_ref_dist
    
    tracks = Tracks(conf=conf)

    print("Initialized ", time.time()-t)
    t = time.time()
    t_all = time.time()

    for read_id,chunk in tracks.iter_reads():

        print(read_id)
        print("    New read ", time.time()-t)
        t = time.time()

        if chunk.any_empty:
            sys.stderr.write(f"Skipping {read_id}\n")
        else:
            #if conf.save:
            #    print(read_id)
            print("    Comparing ", time.time()-t)
            chunk.calc_compare(group_b, calc_jaccard, calc_mean_ref_dist, conf.save)

        print("Read processed ", time.time()-t_all)
        sys.stdout.flush()
        t = time.time()
        t_all = time.time()

DUMP_OPTS = (
    Opt("db_in", "tracks.io", type=str),
    Opt("layers", nargs="+",  help="Layers to retrieve or compute"),
    Opt(("-R", "--ref-bounds"), "tracks", type=ref_coords),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
)

def dump(conf):
    """Output DTW alignment paths and statistics"""


    tracks = Tracks(conf=conf)
    #TODO add layer dependencies (compare requires start/length)
    tracks.set_layers(["start", "length", "bcaln.start", "bcaln.length"] + conf.layers)

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
            sys.stdout.write(track.layers.to_csv(sep="\t", header=False))
        t = time.time()

SUBCMDS = [
    (compare, COMPARE_OPTS),
    (dump, DUMP_OPTS)
]
