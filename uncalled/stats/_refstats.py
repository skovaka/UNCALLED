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
from ..sigproc import ProcRead
from ..argparse import Opt, comma_split
from ..index import BWA_OPTS, str_to_coord
from ..fast5 import Fast5Reader
from ..dtw import TrackIO, LAYER_META

class RefstatsParams(config.ParamGroup):
    _name = "refstats"
RefstatsParams._def_params(
    ("layers", ["current", "dwell"], None, "Layers over which to compute summary statistics"),
    ("stats", ["mean"], None, "Summary statistics to compute for layers specified in \"stats\""),
)

_DESC_FNS = {
    "cov"  : lambda d: d.nobs,
    "mean" : lambda d: d.mean, 
    "stdv" : lambda d: np.power(d.variance,2), 
    "var"  : lambda d: d.variance, 
    "skew" : lambda d: d.skewness, 
    "kurt" : lambda d: d.kurtosis,
    "min"  : lambda d: d.minmax.min, 
    "max"  : lambda d: d.minmax.max}

_AGG_FNS = {
    "mean" : np.mean, 
    "stdv" : np.std, 
    "var"  : np.var,
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : np.min, 
    "max"  : np.min
}

LAYER_STATS = {"min", "max", "mean", "stdv", "var", "skew", "kurt"}
COMPARE_STATS = {"ks"}
ALL_STATS = LAYER_STATS | COMPARE_STATS

class SplitStats:
    def __init__(self, stats, track_count):
        self.layer = [s for s in stats if s in LAYER_STATS]
        self.compare = [s for s in stats if s in COMPARE_STATS]

        self.layer_agg = [_AGG_FNS[s] for s in self.layer]

        if len(self.layer) + len(self.compare) != len(stats):
            bad_stats = [s for s in stats if s not in ALL_STATS]
            raise ValueError("Unknown stats: %s (options: %s)" % (", ".join(bad_stats), ", ".join(ALL_STATS)))

        if len(self.compare) > 0 and track_count != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(self.compare))

class _Refstats:

    def __call__(self, tracks, *args, **kwargs):
        conf, prms = config._init_group("refstats", *args, **kwargs)

        if not isinstance(prms.stats, SplitStats):
            prms.stats = SplitStats(prms.stats, len(tracks))

        stats = dict()
        grouped = [t.layers["dtw"][prms.layers].groupby(level="mref") for t in tracks]

        for track,groups in zip(tracks, grouped):
            stats[track.name] = groups.agg(prms.stats.layer_agg)
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
        stats.index = coords.mref_to_ref_index(stats.index, multi=True)
        #print(stats)
        #print(nt.kmer_to_str(coords.kmers))
        #stats.insert(0, "kmer", nt.kmer_to_str(coords.kmers))

        return stats.dropna()


    @staticmethod
    def describe(track, layers, stats):
        df = {}
        for layer in layers:
            desc = scipy.stats.describe(track[layer], axis=0, nan_policy="omit")
            for stat in stats:
                name = ".".join([layer, stat])
                fn = _Refstats._DESC_FNS[stat]
                df[name] = fn(desc)
                
        return pd.DataFrame(df, index=track.ref_coords.index)
        
    @staticmethod
    def ks(track_a=None, track_b=None):
        grouped = [t.layers["dtw"][prms.layers].groupby(level="mref") for t in (track_a, track_b)]
        cmps = {l : defaultdict(list) for l in prms.layers}
        for (mref,track_a),(_,track_b) in zip(*grouped):
            mrefs.append(mref)
            for layer in prms.layers:
                a = track_a[layer]
                b = track_b[layer]
                for stat in compare_stats:
                    cmps[layer][stat].append(
                        scipy.stats.stats.ks_2samp(a,b,mode="asymp")[0]
                    )
        stats["cmp"] = pd.concat({k : pd.DataFrame(index=mrefs, data=c) for k,c in cmps.items()}, axis=1) 

        df = {}

        for layer in layers:
            layer_stats = np.zeros(track_a.width)
            for i,rf in enumerate(track_a.ref_coords.index):
                a = track_a[layer,:,rf]
                b = track_b[layer,:,rf]
                ks = scipy.stats.stats.ks_2samp(a,b,mode="asymp")
                layer_stats[i] = ks[0]
            df[layer + ".ks"] = layer_stats

        return pd.DataFrame(df, index=track_a.ref_coords.index)

refstats = _Refstats()

OPTS = (
    Opt("layers", "refstats", type=comma_split,
        help="Comma-separated list of layers over which to compute summary statistics {%s}" % ",".join(LAYER_META.keys())),
    Opt("stats", "refstats", type=comma_split,
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_STATS)),
    Opt("input", "track_io", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "track_io", type=str_to_coord),
    Opt(("-C", "--ref-chunksize"), "track_io"),
)

def main(*args, **kwargs):
    """Summarize and compare DTW stats over reference coordinates"""
    conf, prms = config._init_group("refstats", *args, **kwargs)

    t0 = time.time()

    io = TrackIO(conf=conf)
    conf = io.conf

    if not isinstance(prms.stats, SplitStats):
        prms.stats = SplitStats(prms.stats, io.input_count)
    conf.refstats = prms

    columns = ["ref_name", "ref", "strand", "kmer"]
    for track in io.input_tracks:
        name = track.name
        columns.append(".".join([track.name, "cov"]))
        for layer in prms.layers:
            for stat in prms.stats.layer:
                columns.append(".".join([track.name, layer, stat]))

    for layer in prms.layers:
        for stat in prms.stats.compare:
            columns.append(".".join(["cmp", layer, stat]))

    print("\t".join(columns))

    for coords,tracks in io.iter_refs():
        stats = refstats(tracks, conf=conf)
        sys.stdout.write(stats.to_csv(sep="\t",header=False,na_rep=0))
