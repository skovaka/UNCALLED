import sys, os
import numpy as np
import argparse
import re
import time
import types
import pandas as pd
import scipy.stats
from collections import defaultdict

from .. import config, nt
from ..dtw.tracks import RefstatsSplit, ALL_REFSTATS
from ..sigproc import ProcRead
from ..argparse import Opt, comma_split
from ..index import BWA_OPTS, str_to_coord
from ..fast5 import Fast5Reader


OPTS = (
    Opt("input", "tracks", nargs="+", type=str),
    Opt("refstats_layers", "tracks", type=comma_split,
        help="Comma-separated list of layers over which to compute summary statistics"),# {%s}" % ",".join(LAYERS.keys())),
    Opt("refstats", "tracks", type=comma_split,
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_REFSTATS)),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-C", "--ref-chunksize"), "tracks"),
    Opt(("-c", "--cov"), action="store_true", help="Output track coverage for each reference position"),
    Opt(("-v", "--verbose-refs"), action="store_true", help="Output reference name and strand"),
)

def main(conf):
    """Calculate per-reference-coordinate statistics"""
    from ..dtw import Tracks

    t0 = time.time()

    tracks = Tracks(conf=conf)
    conf = tracks.conf

    stats = RefstatsSplit(conf.tracks.refstats, len(tracks.alns))

    if conf.verbose_refs:
        columns = ["ref_name", "ref", "strand"]
    else:
        columns = ["ref"]

    for track in tracks.alns:
        name = track.name
        if conf.cov:
            columns.append(".".join([track.name, "cov"]))
        for group, layer in conf.tracks.refstats_layers:
            for stat in stats.layer:
                columns.append(".".join([track.name, group, layer, stat]))

    for group,layer in conf.tracks.refstats_layers:
        for stat in stats.compare:
            columns.append(".".join([stat, group, layer, "stat"]))

    print("\t".join(columns))

    for chunk in tracks.iter_refs():
        stats = chunk.calc_refstats(conf.verbose_refs, conf.cov)
        if conf.verbose_refs:
            stats.index = pd.MultiIndex.from_product([
                [chunk.coords.ref_name], stats.index, ["+" if chunk.coords.fwd else "-"]
            ])
        sys.stdout.write(stats.to_csv(sep="\t",header=False,na_rep=0))
