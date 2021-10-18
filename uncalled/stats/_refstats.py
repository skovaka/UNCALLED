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
from ..dtw.track_io import RefstatsSplit, ALL_REFSTATS
from ..sigproc import ProcRead
from ..argparse import Opt, comma_split
from ..index import BWA_OPTS, str_to_coord
from ..fast5 import Fast5Reader

class _Refstats:

    def __call__(self, tracks, *args, **kwargs):
        conf, prms = config._init_group("refstats", *args, **kwargs)

        if not isinstance(prms.stats, SplitStats):
            prms.stats = SplitStats(prms.stats, len(tracks))

        stats = dict()
        grouped = [t.layers["dtw"][prms.layers].groupby(level="mref") for t in tracks]

        for track,groups in zip(tracks, grouped):
            stats[track.name] = groups.agg(prms.stats.layer_agg)
            if prms.cov:
                stats[track.name].insert(0, "cov", groups.size())

        if len(prms.stats.compare) > 0:
            groups_a, groups_b = grouped
            mrefs_a = tracks[0].layers.index.unique("mref")
            mrefs_b = tracks[1].layers.index.unique("mref")
            mrefs = mrefs_a.intersection(mrefs_b)
            cmps = {l : defaultdict(list) for l in prms.layers}
            for mref in mrefs:
                track_a = groups_a.get_group(mref)
                track_b = groups_b.get_group(mref)
                for layer in prms.layers:
                    a = track_a[layer]
                    b = track_b[layer]
                    for stat in prms.stats.compare:
                        cmps[layer][stat].append(
                            scipy.stats.stats.ks_2samp(a,b,mode="asymp")[0]
                        )
            stats["cmp"] = pd.concat({k : pd.DataFrame(index=mrefs, data=c) for k,c in cmps.items()}, axis=1) 

        stats = pd.concat(stats, axis=1, names=["track", "layer", "stat"])

        coords = tracks[0].coords
        stats.index = coords.mref_to_ref_index(stats.index, multi=prms.verbose_coords)

        return stats.dropna()

refstats = _Refstats()


OPTS = (
    Opt("refstats_layers", "track_io", type=comma_split,
        help="Comma-separated list of layers over which to compute summary statistics"),# {%s}" % ",".join(LAYERS.keys())),
    Opt("refstats", "track_io", type=comma_split,
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_REFSTATS)),
    Opt("input", "track_io", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "track_io", type=str_to_coord),
    Opt(("-C", "--ref-chunksize"), "track_io"),
    Opt(("-c", "--cov"), action="store_true"),
    Opt(("-v", "--verbose-refs"), action="store_true"),
)

def main(conf):
    """Summarize and compare DTW stats over reference coordinates"""
    from ..dtw import Tracks

    t0 = time.time()

    io = Tracks(conf=conf)
    conf = io.conf

    stats = RefstatsSplit(conf.track_io.refstats, len(io.aln_tracks))

    if conf.verbose_refs:
        columns = ["ref_name", "ref", "strand"]
    else:
        columns = ["ref"]

    for track in io.aln_tracks:
        name = track.name
        if conf.cov:
            columns.append(".".join([track.name, "cov"]))
        for layer in conf.track_io.refstats_layers:
            for stat in stats.layer:
                columns.append(".".join([track.name, layer, stat]))

    for layer in conf.track_io.refstats_layers:
        for stat in stats.compare:
            columns.append(".".join([stat, layer, "stat"]))

    print("\t".join(columns))

    for coords,tracks in io.iter_refs():
        stats = io.calc_refstats(conf.verbose_refs, conf.cov)
        sys.stdout.write(stats.to_csv(sep="\t",header=False,na_rep=0))
