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

class _Dtwstats:
    LAYERS = set(LAYERS.keys())

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("dtwstats", *args, **kwargs)
        
        conf.track.layers = prms.layers

        tracks = _load_tracks(prms.tracks, conf)

        if not isinstance(tracks, list):
            tracks = [tracks]

        if not isinstance(prms.layers, list):
            prms.layers = [prms.layers]
        
        layer_stats = [s for s in prms.layers if s in self.LAYERS]

        if len(layer_stats) != len(prms.layers):
            bad_stats = [s for s in prms.layers if s not in self.LAYERS]
            raise ValueError("Unknown layers: " + ", ".join(bad_stats))

        layers = dict()

        for stat in layer_stats:
            for track in tracks:
                name = ".".join([os.path.basename(track.prms.path), stat])
                fn = getattr(self, stat, None)
                if fn is None:
                    continue
                fn(track)

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


OPTS = (
    Opt("input", "tracks", nargs="+", type=str),
    Opt("layers", nargs="+",  help="Layers to retrieve or compute"),
    Opt(("-R", "--ref-bounds"), "tracks", type=ref_coords),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
)

dtwstats = _Dtwstats()

def main(conf):
    """Output DTW alignment paths, statistics, and comparisons"""

    tracks = Tracks(conf=conf)
    #TODO add layer dependencies (compare requires start/length)
    tracks.set_layers(["start", "length", "bcaln.start", "bcaln.length"])# + conf.layers)

    layer_groups = {group for group,_ in parse_layers(conf.layers)}

    need_cmp = "cmp" in layer_groups
    need_bc_cmp = "bc_cmp" in layer_groups

    for read_id,tracks in tracks.iter_reads():

        if need_cmp and not np.any([t.has_group("cmp") for t in tracks]):
            if len(tracks) != 2:
                raise ValueError("\"cmp\" can only be computed with two alignment tracks")
            if len(tracks[0].alignments) == 0 or len(tracks[1].alignments) == 0:
                continue
            tracks[0].cmp(tracks[1], write=False)

            tracks = [tracks[0]]

        if need_bc_cmp and not np.all([t.has_group("bc_cmp") for t in tracks]):
            if len(tracks) == 1:
                other = None
            elif len(tracks[1].alignments) > 0:
                other = tracks[1]
            else:
                continue
            print("ERE", tracks[0], other, len(tracks))
            tracks[0].bc_cmp(other, write=True)
            if len(tracks[0].alignments) == 0:
                continue
            tracks = [tracks[0]]
        
        for track in tracks:
            print(track.layers.to_csv(sep="\t"))
