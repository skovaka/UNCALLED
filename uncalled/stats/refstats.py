import sys, os
import numpy as np
import argparse
import re
import time
import types
import pandas as pd
import scipy.stats
from collections import defaultdict

from .. import config
from ..dtw.tracks import RefstatsSplit, ALL_REFSTATS
from ..dtw.layers import parse_layers
from ..argparse import Opt, comma_split
from ..index import str_to_coord

def defstats_single(tracks):
    pass

def refstats(conf):
    """Calculate per-reference-coordinate statistics"""
    from ..dtw import Tracks

    t0 = time.time()

    conf.tracks.shared_refs_only = True

    tracks = Tracks(conf=conf)
    conf = tracks.conf
    conf.shared_refs_only = True
    conf.tracks.load_fast5s = False

    #if conf.tracks.io.processes == 1

    stats = RefstatsSplit(conf.refstats, len(tracks.alns))
    layers = list(parse_layers(conf.tracks.layers, False))

    #if conf.verbose_refs:
    columns = ["ref_name", "ref", "strand"]
    #else:
    #    columns = ["ref"]

    for track in tracks.alns:
        name = track.name
        if conf.cov:
            columns.append(".".join([track.name, "cov"]))
        for group, layer in layers:
            for stat in stats.layer:
                columns.append(".".join([track.name, group, layer, stat]))

    for group,layer in layers:
        for stat in stats.compare:
            columns.append(".".join([stat, group, layer, "stat"]))
            columns.append(".".join([stat, group, layer, "pval"]))

    #columns.append("kmer")

    print("\t".join(columns))

    for chunk in tracks.iter_refs():
        chunk.prms.refstats = conf.refstats
        chunk.prms.refstats_layers = layers

        stats = chunk.calc_refstats(conf.cov)
        if stats is None: continue

        stats = pd.concat({chunk.coords.name : stats}, axis=0)\
                  .reset_index(level="seq.fwd")
        stats["seq.strand"] = stats["seq.fwd"].replace({True:"+",False:"-"})
        stats.set_index("seq.strand",append=True,inplace=True)
        del stats["seq.fwd"]
        sys.stdout.write(stats.to_csv(sep="\t",header=False,na_rep=0))
