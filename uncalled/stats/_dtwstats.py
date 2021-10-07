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
from ..dtw import TrackIO, LAYERS

#class DtwstatsParams(config.ParamGroup):
#    _name = "dtwstats"
#DtwstatsParams._def_params(
#    ("layers", ["model_diff"], None, "Which layers to retrieve or compute"),
#    ("tracks", None, None, "DTW alignment tracks"),
#    #("store", True, bool, "Will store layer in track, otherwise just return new layer"),
#)


class _Dtwstats:
    LAYERS = set(LAYERS.keys())
    #COMPARE_STATS = {"cmp"}
    #ALL_STATS = LAYER_STATS | COMPARE_STATS

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("dtwstats", *args, **kwargs)
        
        conf.track.layers = prms.layers

        tracks = _load_tracks(prms.tracks, conf)

        if not isinstance(tracks, list):
            tracks = [tracks]

        if not isinstance(prms.layers, list):
            prms.layers = [prms.layers]
        
        layer_stats = [s for s in prms.layers if s in self.LAYERS]
        #compare_stats = [s for s in prms.layers if s in self.COMPARE_STATS]

        #if len(layer_stats) + len(compare_stats) != len(prms.layers):
        if len(layer_stats) != len(prms.layers):
            bad_stats = [s for s in prms.layers if s not in self.LAYERS]
            raise ValueError("Unknown layers: " + ", ".join(bad_stats))

        #if len(compare_stats) > 0 and len(tracks) != 2:
        #    raise ValueError("\"%s\" layers can only be computed using exactly two tracks" % "\", \"".join(cmp_stats))

        layers = dict()

        for stat in layer_stats:
            for track in tracks:
                name = ".".join([os.path.basename(track.prms.path), stat])
                fn = getattr(self, stat, None)
                if fn is None:
                    continue
                fn(track)

        #for stat in compare_stats:
        #    fn = getattr(self, stat, None)
        #    if fn is None:
        #        continue
        #    mat = fn(*tracks, prms.store)
        #    layers[stat] = mat
                
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
    #    "dwell" : (
    #        lambda self,a: 1000 * a.aln['length'] / self.conf.read_buffer.sample_rate),
    #    "model_diff" : (
    #    "bcerr" : get_bcerr_layer,

OPTS = (
    Opt("input", "track_io", nargs="+", type=str),
    Opt(("-L", "--layers"), "track_io", type=comma_split, 
        help="Comma-separated list of which layers to retrieve or compute {%s}" % ",".join(_Dtwstats.LAYERS)),
    Opt(("-R", "--ref-bounds"), "track_io", type=ref_coords),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
)

dtwstats = _Dtwstats()

def main(conf):
    """Output DTW alignment paths, statistics, and comparisons"""

    io = TrackIO(conf=conf)

    for coords,tracks in io.iter_refs():
        for track in tracks:
            print(track.layers.to_csv(sep="\t"))

    #for track in tracks:
    #    print("#" + track.prms.path)
    #    for i, read_id in enumerate(track.reads["id"]):
    #        print("##" + read_id)
    #        for j in track.ref_coords.index:
    #            for layer in conf.dtwstats.layers:
    #                sys.stdout.write("%.3f\t"%track[layer,i,j])
    #            sys.stdout.write("\n")
