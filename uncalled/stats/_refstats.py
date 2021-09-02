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
from ..fast5 import Fast5Reader
from ..dtw import Track, ref_coords, LAYER_META
from ..dtw.track import _load_tracks

class RefstatsParams(config.ParamGroup):
    _name = "refstats"
RefstatsParams._def_params(
    ("layers", ["current", "dwell"], list, "Layers over which to compute summary statistics"),
    ("stats", ["mean"], list, "Summary statistics to compute for layers specified in \"stats\""),
    ("tracks", None, None, "DTW alignment tracks"),
)


class _Refstats:
    LAYER_STATS = {"cov", "min", "max", "mean", "stdv", "var", "skew", "kurt"}
    COMPARE_STATS = {"ks"}
    ALL_STATS = LAYER_STATS | COMPARE_STATS

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("refstats", *args, **kwargs)
        
        tracks = _load_tracks(prms.tracks, conf)

        layer_stats = [s for s in prms.stats if s in self.LAYER_STATS]
        compare_stats = [s for s in prms.stats if s in self.COMPARE_STATS]

        if len(layer_stats) + len(compare_stats) != len(prms.stats):
            bad_stats = [s for s in prms.stats if s not in _ALL_STATS]
            raise ValueError("Unknown stats: " + ", ".join(bad_stats))

        if len(compare_stats) > 0 and len(tracks) != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(cmp_stats))
            
        
        dfs = list()
        labels = list()

        if len(layer_stats) > 0:
            for track in tracks:
                df = self.describe(track, prms.layers, layer_stats)
                if len(tracks) > 1:
                    df = df.add_prefix(os.path.basename(track.prms.path)+".")
                dfs.append(df)

        for stat in compare_stats:
            fn = getattr(self, stat)
            df = fn(*tracks, prms.layers)
            dfs.append(df)

        df = pd.concat(dfs, axis=1)

        return df

        #for stat in cmp_stats:
        #    for layer in prms.layers:
        #        df = fn(*tracks, layer)
        #        #track_dfs.append(fn(self, *tracks, layers))

    _DESC_FNS = {
        "cov"  : lambda d: d.nobs,
        "mean" : lambda d: d.mean, 
        "stdv" : lambda d: np.power(d.variance,2), 
        "var"  : lambda d: d.variance, 
        "skew" : lambda d: d.skewness, 
        "kurt" : lambda d: d.kurtosis,
        "min"  : lambda d: d.minmax.min, 
        "max"  : lambda d: d.minmax.max}

    @staticmethod
    def describe(track, layers, stats):
        df = {}
        for layer in layers:
            desc = mstats.describe(track[layer], axis=0)
            for stat in stats:
                name = ".".join([layer, stat])
                fn = _Refstats._DESC_FNS[stat]
                df[name] = fn(desc)
                
        return pd.DataFrame(df, index=track.ref_coords.index)
        
    @staticmethod
    def ks(track_a, track_b, layers):
        df = {}

        for layer in layers:
            layer_stats = np.zeros(track_a.width)
            for i,rf in enumerate(track_a.ref_coords.index):
                a = track_a[layer,:,rf]
                b = track_b[layer,:,rf]
                ks = scipy.stats.mstats.ks_2samp(a,b,mode="asymp")
                layer_stats[i] = ks[0]
            df[layer + ".ks"] = layer_stats

        return pd.DataFrame(df, index=track_a.ref_coords.index)

refstats = _Refstats()

OPTS = (
    Opt("layers", "refstats", type=comma_split,
        help="Comma-separated list of layers over which to compute summary statistics {%s}" % ",".join(LAYER_META.keys())),
    Opt("stats", "refstats", type=comma_split,
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(_Refstats.ALL_STATS)),
    Opt("tracks", "refstats", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "track", type=ref_coords, required=True),
)

def main(*args, **argv):
    """Summarize and compare DTW stats over reference coordinates"""
    df = refstats(*args, **argv)
    print(df.to_csv(sep="\t"))
