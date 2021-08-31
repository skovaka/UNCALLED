import sys, os
import numpy as np
import argparse
import re
import time
import scipy.stats.mstats as mstats
import types
import pandas as pd
import scipy.stats
from sklearn.decomposition import PCA

from .. import config
from ..sigproc import ProcRead
from ..argparse import Opt, comma_split, ref_coords
from ..index import BWA_OPTS
from ..fast5 import Fast5Reader
from ..dtw import Track, ref_coords, LAYER_META
from ..dtw.track import _load_tracks

class ReadstatsParams(config.ParamGroup):
    _name = "readstats"
ReadstatsParams._def_params(
    ("stats", ["mean"], None, "Which statistics to compute and output"),
    ("tracks", None, None, "DTW Alignment Track(s)"),
    ("pca_layer", "current", str, "Which statistics to use for PCA"),
    ("pca_components", 2, int, "Number of principle components to output for the \"pca\" command."),
    ("summary_stats", ["mean"], None, "Summary statistics to compute for \"model_diff\" command."),
)

OPTS = (
    Opt("stats", "readstats", type=comma_split),
    Opt("tracks", "readstats", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "track", type=ref_coords, required=True),
    Opt(("-p", "--pca-components"), "readstats"),
    Opt(("-L", "--pca-layer"), "readstats"),
    Opt(("-s", "--summary-stats"), "readstats", type=comma_split),
)

class _Readstats:
    STATS = {"model_diff", "pca",} #"speed", "hierclust", "kmeans"

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("readstats", *args, **kwargs)
        
        tracks = _load_tracks(prms.tracks, conf, True)

        if isinstance(prms.stats, list):
            stats = prms.stats
        else:
            stats = [prms.stats]

        dfs = list()
        
        for track in tracks:
            for stat in stats:
                fn = getattr(self, stat)
                df = fn(track, prms=prms)
                if len(tracks) > 1:
                    df = df.add_prefix(os.path.basename(track.prms.path)+".")
                dfs.append(df)

        return pd.concat(dfs, axis=1)

    _DESC_FNS = {
        "length" : lambda d: d.mean, 
        "mean" : lambda d: d.mean, 
        "stdv" : lambda d: np.power(d.variance,2), 
        "var"  : lambda d: d.variance, 
        "skew" : lambda d: d.skewness, 
        "kurt" : lambda d: d.kurtosis}

    #@staticmethod
    #def _summary_stats(track, layer, stats)

    @staticmethod
    def model_diff(track, summary_stats=None, prms=None):
        if summary_stats is None:
            if prms is None:
                raise ValueError("Must specify summary_stats or ReadstatsParams for readstats.model_diff")
            summary_stats = prms.summary_stats

        layer = "model_diff"

        desc = mstats.describe(track[layer], axis=1)

        df = pd.DataFrame(index=track.reads["id"])
        for stat in summary_stast:
            name = ".".join(layer, stat)
            df[name] = _DESC_FNS[stat](desc)
            
        return df

    @staticmethod
    def pca(track, layer=None, components=True, prms=None):
        if components is None or layer is None:
            if prms is None:
                raise ValueError("Insufficent parameters for readstats.pca")
            layer = prms.pca_layer
            components = prms.pca_components

        x = track[layer,:,:].T
        pc = PCA(n_components=components).fit_transform(x)
        df = pd.DataFrame(index=track.reads["id"])
        for c in range(components):
            print(pc[c,:])
            df["pc%d" % (c+1)] = pc[:,c]

        return df

readstats = _Readstats()

def main(*args, **argv):
    """Outputs statistics for each reference position over one or more tracks"""
    df = readstats(*args, **argv)
    print(df.to_csv(sep="\t"))
