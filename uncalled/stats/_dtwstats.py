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
from ..dtw import Track, ref_coords, LAYER_META
from ..dtw.track import _load_tracks

class DtwstatsParams(config.ParamGroup):
    _name = "dtwstats"
DtwstatsParams._def_params(
    ("stats", ["model_diff"], list, "Which statistics to compute as new track layers"),
    ("tracks", None, None, "DTW Alignment Track(s)"),
    ("store", True, bool, "Will store layer in track, otherwise just return new layer"),
)

OPTS = (
    Opt("stats", "dtwstats", type=comma_split),
    Opt("tracks", "dtwstats", nargs="+", type=str),
    Opt(("-R", "--ref-bounds"), "track", type=ref_coords, required=True),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
)

class _Dtwstats:
    LAYER_STATS = set(LAYER_META.keys())
    COMPARE_STATS = {"overlap"}
    ALL_STATS = LAYER_STATS | COMPARE_STATS

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("dtwstats", *args, **kwargs)
        
        conf.track.layers = prms.stats

        tracks = _load_tracks(prms.tracks, conf)

        if not isinstance(tracks, list):
            tracks = [tracks]

        if not isinstance(prms.stats, list):
            prms.stats = [prms.stats]
        
        layer_stats = [s for s in prms.stats if s in self.LAYER_STATS]
        compare_stats = [s for s in prms.stats if s in self.COMPARE_STATS]

        if len(layer_stats) + len(compare_stats) != len(prms.stats):
            bad_stats = [s for s in prms.stats if s not in self.ALL_STATS]
            raise ValueError("Unknown stats: " + ", ".join(bad_stats))

        if len(compare_stats) > 0 and len(tracks) != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(cmp_stats))

        layers = dict()

        for stat in layer_stats:
            for track in tracks:
                name = ".".join([os.path.basename(track.prms.path), stat])
                fn = getattr(self, stat, None)
                if fn is None:
                    continue
                mat = fn(track, prms.store)
                layers[name] = mat

        for stat in compare_stats:
            fn = getattr(self, stat, None)
            if fn is None:
                continue
            mat = fn(*tracks, prms.store)
            layers[stat] = mat
                
        return tracks #layers

    @staticmethod
    def _add_layer(track, name, layer, store):
        if store:
            track.add_layer(name, layer)
            return track[name]
        return layer

    @staticmethod
    def model_diff(track, store=True):
        layer = track["current"] - track.model[track["kmer"].astype(int)]
        return _Dtwstats._add_layer(track, "model_diff", layer, store)

    @staticmethod
    def dwell(track, store=True):
        layer = 1000 * track['length'] / track.conf.read_buffer.sample_rate
        return _Dtwstats._add_layer(track, "dwell", layer, store)
        
    #@staticmethod
    #def bcerr(self, aln):
    #    bcerr = aln.bcerr#.reindex(aln.aln.index)
    #    ret = pd.Series(np.nan, bcerr.index)
    #    subs = bcerr[bcerr["type"]=="SUB"]
    #    ret[subs.index] = subs["seq"].replace({"A":0,"C":1,"G":2,"T":3})
    #    ret[bcerr["type"]=="DEL"] = 4
    #    ret[bcerr["type"]=="INS"] = 5
    #    return ret

    #LAYER_FNS = {
    #    #"id" : (
    #    #    lambda self,a: a.id),
    #    "kmer" : (
    #        lambda self,a: self.load_aln_kmers(a)),
    #    "current" : (
    #        lambda self,a: a.aln["current"]),
    #    "dwell" : (
    #        lambda self,a: 1000 * a.aln['length'] / self.conf.read_buffer.sample_rate),
    #    "model_diff" : (
    #    "bcerr" : get_bcerr_layer,

dtwstats = _Dtwstats()

def main(conf):
    """Outputs statistics for each reference position over one or more tracks"""
    tracks = dtwstats(conf=conf)

    layers = conf.dtwstats.stats

    for track in tracks:
        print("#" + track.prms.path)
        for i, read_id in enumerate(track.reads["id"]):
            print("##" + read_id)
            for j in track.ref_coords.index:
                for layer in layers:
                    sys.stdout.write("%.3f\t"%track[layer,i,j])
                sys.stdout.write("\n")
