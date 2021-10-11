#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict, namedtuple
import re
import time
from typing import NamedTuple
from matplotlib.colors import Normalize
import pandas as pd
import copy
import sqlite3

import scipy.stats
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt

from ..pafstats import parse_paf, PafEntry
from ..argparse import Opt, ref_coords
from .. import nt, PoreModel, config, index 
from ..index import load_index, RefCoord

KMER_LAYER       = "kmer"
CURRENT_LAYER    = "current"
DWELL_LAYER      = "dwell"
MODEL_DIFF_LAYER = "model_diff"
DEFAULT_LAYERS = [CURRENT_LAYER, DWELL_LAYER, MODEL_DIFF_LAYER]

LayerMeta = namedtuple("LayerMeta", ["type", "label", "fn"])

#TODO probably move this to AlnTrack
LAYERS = {
    "ref" : {
        "coord" : LayerMeta(int, "Reference Coordinate", 
            lambda track: track.coords.mref_to_ref(track.layer_mrefs)),
        "name" : LayerMeta(str, "Reference Name",
            lambda track: [track.coords.ref_name]*len(track.layers)),
        "fwd" : LayerMeta(bool, "Is on fwd strand",
            lambda track: [track.coords.fwd]*len(track.layers)),
    }, "dtw" : {
        "start" : LayerMeta(int, "Sample Start", None),
        "length" : LayerMeta(int, "Sample Length", None),
        "current" : LayerMeta(float, "Current (pA)", None),
        "end" : LayerMeta(int, "Sample End",  
            lambda track: track.layers["dtw","start"] + track.layers["dtw","length"]),
        "middle" : LayerMeta(float, "Sample Middle",  
            lambda track: track.layers["dtw","start"] + (track.layers["dtw","length"] / 2)),
        "dwell" : LayerMeta(float, "Dwell Time (ms/nt)",
            lambda track: 1000 * track.layers["dtw","length"] / track.conf.read_buffer.sample_rate),
        "model_diff" : LayerMeta(float, "Model pA Difference",
            lambda track: track.layers["dtw","current"] - track.model[track.kmers]),
        "kmer" : LayerMeta(str, "Reference k-mer",
            lambda track: track.coords.kmers[track.layer_mrefs].to_numpy()),
        "base" : LayerMeta(str, "Reference base",
            lambda track: nt.kmer_base(track.coords.kmers[track.layer_mrefs], 2)),
        "bp" : LayerMeta(int, "Basecaller Base Index", None),
        "err" : LayerMeta(str, "Basecalled Alignment Error", None)}
}

class AlnTrack:

    #TODO get rid of htis
    CMP_LAYERS = ["current", "dwell"]

    def __init__(self, db, track_id, name, desc, groups, conf):
        self.db = db
        #self.index = index
        self.id = track_id
        self.name = name
        self.desc = desc
        self.groups = groups.split(",")
        self.conf = conf

        self.read_ids = set()

        self.mat = None
        self.coords = None
        self.layers = None

        self.model = PoreModel(self.conf.pore_model) 

    def _group_layers(self, group, layers):
        return pd.concat({group : layers}, names=["group", "layer"], axis=1)
        
    def add_layer_group(self, group, df, aln_id=None):
        aln_id = self._default_id(aln_id)
        df = self._group_layers(group, df)#pd.concat({group : layers}, names=["group", "layer"], axis=1)

        df.index = pd.MultiIndex.from_product(
                        [df.index, [aln_id]], 
                        names=["mref", "aln_id"])

        self.db.write_layers(df)

        if self.layers is None:
            self.layers = df
        else:
            self.layers = pd.concat([self.layers, df], axis=1)

    def aln_ref_coord(self, aln_id):
        return RefCoord(*self.alignments[["ref_name","ref_start","ref_end","fwd"]].loc[aln_id])

    def _default_id(self, aln_id):
        if aln_id is None:
            if len(self.alignments) == 1:
                return self.alignments.index[0]
            raise ValueError("Must specify aln_id for Track with more than one alignment loaded")
        return aln_id

    def set_data(self, coords, alignments, layers, load_mat=False):
        self.coords = coords
        self.layers = layers#pd.concat(layers, names=["group", "layer"], axis=1)

        self.alignments = alignments
        self.alignments.sort_values(["fwd", "ref_start"], inplace=True)

        if self.coords.kmers is not None:
            mrefs = self.layers.index.get_level_values("mref")
            if self.coords.stranded:
                self.kmers = self.coords.kmers[mrefs]
            else:
                self.kmers = pd.concat(self.coords.kmers)[mrefs]

        self.has_fwd = np.any(self.alignments['fwd'])
        self.has_rev = not np.all(self.alignments['fwd'])

    def calc_layers(self, layers):
        for group, layer in layers:
            if not group in self.layers or not layer in self.layers[group].columns:
                vals = LAYERS[group][layer].fn(self)
                self.layers[group,layer] = vals

    def load_mat(self):
        df = self.layers.copy()
        df["aln_id"] = df.index.get_level_values("aln_id")
        df.index = self.coords.mref_to_ref_index(df.index.get_level_values("mref"), multi=False)
        df = df.reset_index()

        self.mat = df.pivot(index="aln_id", columns=["ref"]) \
                     .rename_axis(("group","layer","ref"), axis=1) \
                     .sort_index(axis=1, level="ref")

        self.mat = self.mat.reindex(self.alignments.index, copy=False)

        self.width = len(self.coords.refs)
        self.height = len(self.alignments)

    def sort_coord(self, layer, ref):
        order = (-self.mat[layer,ref].fillna(0)).argsort()
        self.sort(order)

    def sort(self, order):
        self.mat = self.mat.iloc[order]
        self.alignments = self.alignments.iloc[order]

    def get_aln_layers(self, aln_id, group=None, layers=None):
        #self.calc_layers([("re)])
        if layers is not None:
            self.calc_layers([(group, layer) for layer in layers])
        df = self.layers.xs(aln_id, level="aln_id").set_index(self.layer_refs)
        if group is None:
            return df
        if layers is None:
            return df[group]
        return df[group][layers]

    def compare(self, other):
        groups_b = other.alignments.groupby("read_id")

        for id_a, aln_a in self.alignments.iterrows():
            read_id = aln_a["read_id"]
            for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                dtw_a = self.get_aln_layers(id_a, "dtw", ["start","end"])
                dtw_b = other.get_aln_layers(id_b, "dtw", ["start","end"])
                merge = dtw_a.join(dtw_b, on="ref", lsuffix="_a", rsuffix="_b")

                #intv_a = pd.IntervalIndex.from_arrays(dtw_a["start"], dtw_a["end"])
                #intv_b = pd.IntervalIndex.from_arrays(dtw_b["start"], dtw_b["end"])

                starts = merge[["start_a", "start_b"]]
                ends = merge[["end_a", "end_b"]]
                jac = 1 - (
                    (ends.min(axis=1) - starts.max(axis=1)).clip(0) /
                    (ends.max(axis=1) - starts.min(axis=1)))

                self.layers["dtw","jac_dist"] = jac[self.layer_refs].to_numpy()


    @property
    def layer_mrefs(self):
        return self.layers.index.get_level_values("mref")

    @property
    def layer_refs(self):
        return self.coords.mref_to_ref_index(self.layer_mrefs)
