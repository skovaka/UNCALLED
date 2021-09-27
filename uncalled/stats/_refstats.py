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
    ("layers", ["current", "dwell"], list, "Layers over which to compute summary statistics"),
    ("stats", ["mean"], list, "Summary statistics to compute for layers specified in \"stats\""),
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
#    "cov"  : len,
    "mean" : np.mean, 
    "stdv" : np.std, 
    "var"  : np.var,
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : np.min, 
    "max"  : np.min
}

class _Refstats:
    LAYER_STATS = {"min", "max", "mean", "stdv", "var", "skew", "kurt"}
    COMPARE_STATS = {"ks"}
    ALL_STATS = LAYER_STATS | COMPARE_STATS

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("refstats", *args, **kwargs)

        io = TrackIO(conf=conf)
        conf = io.conf

        #tracks = _load_tracks(prms.tracks, conf)

        layer_stats = [s for s in prms.stats if s in self.LAYER_STATS]
        compare_stats = [s for s in prms.stats if s in self.COMPARE_STATS]

        if len(layer_stats) + len(compare_stats) != len(prms.stats):
            bad_stats = [s for s in prms.stats if s not in self.ALL_STATS]
            raise ValueError("Unknown stats: %s (options: %s)" % (", ".join(bad_stats), ", ".join(self.ALL_STATS)))

        if len(compare_stats) > 0 and io.input_count != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(cmp_stats))
            
        #for coords,tracks in io.iter_refs():
        #    for t in tracks:
        #        grouped = t.layers["dtw"][prms.layers].groupby("mref")
        #        desc = grouped.aggregate(scipy.stats.describe)


        layer_agg = [_AGG_FNS[s] for s in layer_stats]

        columns = ["ref_name", "ref", "strand", "kmer"]
        for track in io.input_tracks:
            name = track.name
            columns.append(".".join([track.name, "cov"]))
            for layer in prms.layers:
                for stat in layer_stats:
                    columns.append(".".join([track.name, layer, stat]))
        for layer in prms.layers:
            for stat in compare_stats:
                columns.append(".".join(["cmp", layer, stat]))

        print("\t".join(columns))

        t = time.time()
        for coords,tracks in io.iter_refs():
            stats = dict()
            grouped = [t.layers["dtw"][prms.layers].groupby(level="mref") for t in tracks]
            for track,groups in zip(tracks, grouped):
                stats[track.name] = groups.agg(layer_agg)
                stats[track.name].insert(0, "cov", groups.size())

            if len(compare_stats) > 0:
                mrefs = list()
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


            #for i,rf in enumerate(track_a.ref_coords.index):
            #    a = track_a[layer,:,rf]
            #    b = track_b[layer,:,rf]
            #    ks = scipy.stats.stats.ks_2samp(a,b,mode="asymp")
            #    layer_stats[i] = ks[0]
            #df[layer + ".ks"] = layer_stats
                

            stats = pd.concat(stats, axis=1, names=["track", "layer", "stat"])

            stats.index = coords.mref_to_ref_index(stats.index)
            #pd.MultiIndex.from_product([[coords.ref_name], coords.mref_to_ref(stats.index), ["+" if coords.fwd else "-"]], names=["ref_name", "ref", "strand"])

            stats.insert(0, "kmer", nt.kmer_to_str(coords.kmers))
            sys.stdout.write(stats.to_csv(sep="\t",header=False))

        print("#",time.time()-t)

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
    def ks(track_a, track_b, layers):
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
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(_Refstats.ALL_STATS)),
    Opt("input", "track_io", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "track_io", type=str_to_coord),
    Opt(("-C", "--ref-chunksize"), "track_io"),
)

def main(*args, **argv):
    """Summarize and compare DTW stats over reference coordinates"""
    df = refstats(*args, **argv)
    #print(df.to_csv(sep="\t"))
