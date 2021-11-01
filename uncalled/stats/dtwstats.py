"""Compute new DTW alignemnt layers

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
    Opt("input", "tracks", nargs="+", type=str),
    Opt(("-b", "--bcaln"), action="store_true", help="Compare against basecalled alignment. If two tracks input will look for \"bcaln\" group in second track, otherwise will look in the first track."),
    Opt(("-s", "--save"), action="store_true", help="Will save in database if included, otherwise will output to TSV"),
    #Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
)

def compare(conf):
    group_b = "bcaln" if conf.bcaln else "dtw"

    if conf.bcaln:
        conf.tracks.layers += [("bc_cmp", "mean_ref_dist")] #LAYERS["bc_cmp"]["mean_ref_dist"].deps
    else:
        conf.tracks.layers += [("cmp", "mean_ref_dist")] #LAYERS["cmp"]["mean_ref_dist"].deps
    
    tracks = Tracks(conf=conf)
    for read_id,_ in tracks.iter_reads():
        if tracks.any_empty:
            sys.stderr.write(f"Skipping {read_id}\n")
        else:
            if conf.save:
                print(read_id)
            tracks.calc_compare(group_b, conf.save)

SUBCMDS = [
    (compare, COMPARE_OPTS)
]

#OPTS = (
#    Opt("input", "tracks", nargs="+", type=str),
#    Opt("layers", nargs="+",  help="Layers to retrieve or compute"),
#    Opt(("-R", "--ref-bounds"), "tracks", type=ref_coords),
#    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
#    Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
#)
#
#dtwstats = _Dtwstats()
#
#def main(conf):
#    """Output DTW alignment paths, statistics, and comparisons"""
#
#    tracks = Tracks(conf=conf)
#    #TODO add layer dependencies (compare requires start/length)
#    tracks.set_layers(["start", "length", "bcaln.start", "bcaln.length"])# + conf.layers)
#
#    layer_groups = {group for group,_ in parse_layers(conf.layers)}
#
#    need_cmp = "cmp" in layer_groups
#    need_bc_cmp = "bc_cmp" in layer_groups
#
#    for read_id,tracks in tracks.iter_reads():
#
#        if need_cmp and not np.any([t.has_group("cmp") for t in tracks]):
#            if len(tracks) != 2:
#                raise ValueError("\"cmp\" can only be computed with two alignment tracks")
#            if len(tracks[0].alignments) == 0 or len(tracks[1].alignments) == 0:
#                continue
#            tracks[0].cmp(tracks[1], write=False)
#
#            tracks = [tracks[0]]
#
#        if need_bc_cmp and not np.all([t.has_group("bc_cmp") for t in tracks]):
#            if len(tracks) == 1:
#                other = None
#            elif len(tracks[1].alignments) > 0:
#                other = tracks[1]
#            else:
#                continue
#            print("ERE", tracks[0], other, len(tracks))
#            tracks[0].bc_cmp(other, write=True)
#            if len(tracks[0].alignments) == 0:
#                continue
#            tracks = [tracks[0]]
#        
#        for track in tracks:
#            print(track.layers.to_csv(sep="\t"))
