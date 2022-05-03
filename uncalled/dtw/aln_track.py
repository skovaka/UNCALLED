#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict, namedtuple
import re
import time
from typing import NamedTuple
import pandas as pd
import copy

from dataclasses import dataclass
from typing import Callable

import scipy.stats
from sklearn.decomposition import PCA

from ..pafstats import parse_paf, PafEntry
from ..argparse import Opt, ref_coords
from .. import PoreModel, config, index 
from ..index import load_index, RefCoord
from .bcaln import Bcaln

KMER_LAYER       = "kmer"
CURRENT_LAYER    = "current"
DWELL_LAYER      = "dwell"
MODEL_DIFF_LAYER = "model_diff"
DEFAULT_LAYERS = [CURRENT_LAYER, DWELL_LAYER, MODEL_DIFF_LAYER]


LayerMeta = namedtuple("LayerMeta", ["type", "label", "fn", "deps"], defaults=[None,None])

#@dataclass 
#class LayerMeta:
#    _type : type
#    label : str
#    fn    : Callable = None
#    deps  : list = None


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
        "start" : LayerMeta(int, "Sample Start"),
        "length" : LayerMeta(int, "Sample Length"),
        "current" : LayerMeta(float, "Current (pA)"),
        "end" : LayerMeta(int, "Sample End",  
            lambda track: track.layers["dtw","start"] + track.layers["dtw","length"],
            [("dtw", "start"), ("dtw", "length")]),
        "middle" : LayerMeta(float, "Sample Middle",  
            lambda track: track.layers["dtw","start"] + (track.layers["dtw","length"] / 2),
            [("dtw", "start"), ("dtw", "length")]),
        "dwell" : LayerMeta(float, "Dwell (ms/nt)",
            lambda track: 1000 * track.layers["dtw","length"] / track.conf.read_buffer.sample_rate,
            [("dtw", "length")]),
        "model" : LayerMeta(float, "Model Current",
            lambda track: track.model[track.kmers],
            [("dtw", "current")]),
        "model_diff" : LayerMeta(float, "Model pA Diff.",
            lambda track: track.layers["dtw","current"] - track.model[track.kmers],
            [("dtw", "current")]),
        "kmer" : LayerMeta(str, "Reference k-mer",
            lambda track: track.kmers),
        "base" : LayerMeta(str, "Reference base",
            lambda track: track.model.kmer_base(track.kmers, 2)),
    }, "bcaln" : {
        "start" : LayerMeta(int, "BC Sample Start"),
        "length" : LayerMeta(int, "BC Sample Length"),
        "end" : LayerMeta(int, "Sample End",  
            lambda track: track.layers["bcaln","start"] + track.layers["bcaln","length"],
            [("bcaln", "start"), ("bcaln", "length")]),
        "middle" : LayerMeta(float, "Sample Middle",  
            lambda track: track.layers["bcaln","start"] + (track.layers["bcaln","length"] / 2),
            [("bcaln", "start"), ("bcaln", "length")]),
        "bp" : LayerMeta(int, "Basecaller Base Index"),
        "error" : LayerMeta(str, "Basecalled Alignment Error"),
    }, "cmp" : {
        "jaccard" : LayerMeta(int, "Jaccard Distance", None, 
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length")]),
        "mean_ref_dist" : LayerMeta(int, "Mean Ref. Distance", None,
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length")]),
    }, "bc_cmp" : {
        "jaccard" : LayerMeta(int, "Jaccard Distance", None, 
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length"), ("bcaln", "start"), ("bcaln", "end"), ("bcaln", "length")]),
        "mean_ref_dist" : LayerMeta(int, "Mean Ref. Distance", None,                   
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length"), ("bcaln", "start"), ("bcaln", "end"), ("bcaln", "length")]),
    }
}

def parse_layer(layer):
    
    if isinstance(layer, str):
        spl = layer.split(".")
    elif isinstance(layer, tuple):
        spl = layer
    else:
        raise ValueError("Layer must be string or tuple")

    if len(spl) == 2:
        group,layer = spl
    #TODO allow for full group specification
    elif len(spl) == 1:
        if layer in LAYERS:
            group = layer
            for layer in LAYERS[group].keys():
                yield (group, layer)
            return
        group = "dtw"
        layer = spl[0]
    else:
        raise ValueError("Invalid layer specifier \"{layer}\", must contain at most one \".\"")

    if not group in LAYERS:
        opts = ",".join(LAYERS.keys())
        raise ValueError("Invalid layer group \"{group}\". Options: {opts}")

    group_layers = LAYERS[group]

    if not layer in group_layers:
        opts = "\", \"".join(group_layers.keys())
        raise ValueError(f"Invalid layer \"{group}.{layer}\". Options: \"{opts}\"")

    yield (group, layer)

def parse_layers(layers, add_deps=True):
    if layers is None:
        return None

    db_layers = list() 
    fn_layers = list() 

    if isinstance(layers, str):
        layers = layers.split(",")

    ret = list()

    parsed = set()

    for layerstr in layers:
        for layer in parse_layer(layerstr):
            if not layer in parsed:
                parsed.add(layer)

                if add_deps:
                    deps = LAYERS[layer[0]][layer[1]].deps
                    if deps is not None:
                        for dep in deps:
                            if not dep in parsed:
                                parsed.add(dep)
                                yield dep

                yield layer

def add_group_level(self, group, layers):
    return pd.concat({group : layers}, names=["group", "layer"], axis=1)

class AlnTrack:
    def __init__(self, *args, **kwargs):
        if isinstance(args[0], AlnTrack):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

        self.mat = None
        #self.coords = None
        #self.alignments = pd.DataFrame()
        #self.layers = pd.DataFrame()

    def _init_new(self, db, track_id, name, desc, conf, model, fast5s=None):
        self.db = db
        self.id = track_id
        self.name = name
        self.desc = desc
        self.conf = conf

        self.fast5s = fast5s #TODO get rid of this
        self.model = model 

    def _init_slice(self, p, coords=None, alignments=None, layers=None, order=["fwd", "ref_start"]):
        self._init_new(p.db, p.id, p.name, p.desc, p.conf, p.model, p.fast5s)
        self.set_data(coords, alignments, layers, order)

    def set_data(self, coords, alignments, layers, order=["fwd", "ref_start"]):

        self.coords = coords
        self.alignments = alignments
        self.layers = layers

        if not self.coords.stranded and (self.all_fwd or self.all_rev):
            self.coords = self.coords.ref_slice(fwd=self.all_fwd)

        isnone = [coords is None, alignments is None, layers is None]
        if np.all(isnone) or len(alignments) == 0:
            return
        elif np.any(isnone):
            raise ValueError("Must specify AlnTrack coords, alignments, and layers")
        print(self.alignments)
        print(order)
        self.alignments = self.alignments.sort_values(order)

        if self.layers.index.names[0] == "mref":
            self.layers = self.layers.rename(index=coords.mref_to_ref, level=0)
            self.layers.index.names = ("ref", "aln_id")

        self.layers = self.layers.sort_index()

        self.layer_fwds = self.alignments.loc[self.layer_aln_ids, "fwd"].to_numpy()

        if self.coords.ref_kmers is not None:
            kidx = pd.MultiIndex.from_arrays([self.layer_fwds, self.layer_refs])
            self.kmers = self.coords.ref_kmers.reindex(kidx)
            self.kmers.index = self.layers.index

        self.has_fwd = np.any(self.alignments['fwd'])
        self.has_rev = not np.all(self.alignments['fwd'])

    @property
    def all_fwd(self):
        return np.all(self.alignments["fwd"])

    @property
    def all_rev(self):
        return not np.any(self.alignments["fwd"])

    #def slice(self, ref_start=0, ref_end=np.inf, aln_ids=None):
    def slice(self, coords=slice(None), aln_ids=None, reads=None, order=["fwd","ref_start"]):
        #ref_start = max(self.layer_refs.min(), ref_start)
        #ref_end = min(self.layer_refs.max()+1, ref_end)
        #coords = self.coords.ref_slice(ref_start, ref_end)
        if self.empty:
            return AlnTrack(self, coords, self.alignments, self.layers)
        layer_refs = self.layers.index.get_level_values("ref")

        layers = self.layers.loc[layer_refs.isin(coords.refs)]

        if reads is not None:
            if aln_ids is not None:
                raise ValueError("Only one of 'aln_ids' and 'reads' can be specified")
            aln_ids = self.alignments.index[self.alignments["read_id"].isin(reads)]

        if aln_ids is not None: 
            layer_alns = layers.index.get_level_values("aln_id")
            aln_ids = layer_alns.unique().intersection(aln_ids)
            layers = layers.loc[layer_alns.isin(aln_ids)] 
        else:
            aln_ids = slice(None)

        alignments = self.alignments.loc[aln_ids]

        return AlnTrack(self, coords, alignments, layers, 
                        order=order) #TODO should handle order in tracks
        

    @property
    def empty(self):
        return self.coords is None or len(self.alignments) == 0

    def _group_layers(self, group, layers):
        return pd.concat({group : layers}, names=["group", "layer"], axis=1)

    def _aln_id_or_default(self, aln_id):
        if aln_id is None:
            if len(self.alignments) == 1:
                return self.alignments.index[0]
            raise ValueError("Must specify aln_id for Track with more than one alignment loaded")
        return aln_id
        
    def add_layer_group(self, group, layers, aln_id, overwrite):
        aln_id = self._aln_id_or_default(aln_id)
        df = pd.concat({group : layers}, names=["group", "layer"], axis=1)

        df.index = pd.MultiIndex.from_product(
                        [df.index, [aln_id]], 
                        names=["mref", "aln_id"])

        if self.layers is None or overwrite:
            self.layers = df
        else:
            self.layers = pd.concat([self.layers, df], axis=1)

        return df

    def aln_ref_coord(self, aln_id):
        return RefCoord(*self.alignments[["ref_name","ref_start","ref_end","fwd"]].loc[aln_id])

    def has_group(self, group):
        return group in self.layers.columns.get_level_values(0)


    def calc_layers(self, layers):
        for group, layer in layers:
            if not group in self.layers or not layer in self.layers[group].columns:
                meta = LAYERS[group][layer]

                #Make sure layer dependencies exist
                if not self.empty and (meta.deps is None or len(self.layers.columns.intersection(meta.deps)) == len(meta.deps)):
                    fn = meta.fn
                    if fn is None:
                        raise ValueError("Layer not found: {group}.{layer}")
                    vals = fn(self)
                    self.layers[group,layer] = vals
                else:
                    self.layers[group,layer] = pd.NA

    def load_mat(self):
        df = self.layers.copy()
        #df["aln_id"] = df.index.get_level_values("aln_id")
        df = df.reset_index()

        self.mat = df.pivot(index="aln_id", columns=["ref"]) \
                     .rename_axis(("group","layer","ref"), axis=1) \
                     .sort_index(axis=1)

        self.mat = self.mat.reindex(self.alignments.index, copy=False)

        self.width = len(self.coords.refs)
        self.height = len(self.alignments)

    def sort_coord(self, layer, ref):
        order = (-self.mat[layer,ref].fillna(0)).argsort()
        self.sort(order)

    def sort(self, order):
        self.mat = self.mat.iloc[order]
        self.alignments = self.alignments.iloc[order]

    def get_aln_layers(self, aln_id, group=None, layers=None, drop_level=True):
        if aln_id not in self.layers.index.get_level_values("aln_id"):
            return None

        if layers is not None:
            self.calc_layers([(group, layer) for layer in layers])

        df = self.layers.xs(aln_id, level="aln_id", drop_level=drop_level)#.set_index(self.layer_refs)

        if group is None:
            return df
        if layers is None:
            return df[group]
        return df[group][layers]

    def cmp(self, other, calc_jaccard, calc_mean_ref_dist):
        groups_b = other.alignments.groupby("read_id")

        df = pd.DataFrame(
            columns=["aln_b", "group_b", "mean_ref_dist", "jaccard"],
            index = self.layers.index
        )

        for id_a, aln_a in self.alignments.iterrows():
            read_id = aln_a["read_id"]
            dtw_a = self.get_aln_layers(id_a, "dtw", ["start","end"], False)

            for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                self._compare_alns(dtw_a, other, id_b, "dtw", df, calc_jaccard, calc_mean_ref_dist)

        df["group_b"] = "dtw"
        return df.set_index(["aln_b", "group_b"], append=True)

    def bc_cmp(self, other=None, write=False):
        if other is not None:
            groups_b = other.alignments.groupby("read_id")

        df = pd.DataFrame(
            columns=["aln_b", "group_b", "mean_ref_dist", "jaccard"],
            index = self.layers.index
        )

        for id_a, aln_a in self.alignments.iterrows():
            dtw = self.get_aln_layers(id_a, "dtw", ["start","end"], False)
            if dtw is None:
                continue

            if other is None:
                self._compare_alns(dtw, self, id_a, "bcaln", df)
            else:
                read_id = aln_a["read_id"]
                for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                    self._compare_alns(dtw, other, id_b, "bcaln", df)

        df["group_b"] = "bcaln"
        return df.set_index(["aln_b", "group_b"], append=True)


    def _compare_alns(self, aln_a, other, id_b, group, df, calc_jaccard=True, calc_mean_ref_dist=True):
        aln_b = other.get_aln_layers(id_b, group, ["start","end"], False) \
                     .reset_index(level="aln_id") \
                     .rename(columns={"aln_id" : "aln_b"})
        #if group == "bcaln":
        #    aln_b.index = aln_b.index+1

        merge = aln_a.join(aln_b, on="ref", lsuffix="_a", rsuffix="_b") \
                     .dropna(how="all")

        #print(merge)
                     
        starts = merge[["start_a", "start_b"]]
        ends = merge[["end_a", "end_b"]]

        df.loc[merge.index, "aln_b"] = merge["aln_b"].astype("Int32")

        if calc_jaccard:
            df.loc[merge.index, "jaccard"] = 1 - (
                (ends.min(axis=1) - starts.max(axis=1)).clip(0) /
                (ends.max(axis=1) - starts.min(axis=1)))
        #else:
        #    jaccard = pd.NA

        if calc_mean_ref_dist:
            intvs_a = pd.IntervalIndex.from_arrays(merge["start_a"], merge["end_a"])
            intvs_b = pd.IntervalIndex.from_arrays(merge["start_b"], merge["end_b"])

            mean_ref_dist = pd.Series(index=merge.index)
            refs = merge.index.get_level_values(0)

            def weighted_dist(idxs, idx_b, swap):
                if swap:
                    a = "b"
                    b = "a"
                else:
                    a = "a"
                    b = "b"
                coord_a = merge.loc[idxs]
                coord_b = merge.loc[idx_b]
                samp_overlaps = coord_a["end_"+a].clip(upper=coord_b["end_"+b]) - coord_a["start_"+a].clip(lower=coord_b["start_"+b])
                ref_dists = np.abs(idxs.get_level_values(0).to_numpy()-idx_b[0])
                return np.sum(samp_overlaps), np.sum(ref_dists * samp_overlaps)

            for i,idx in enumerate(merge.index):
                ref,_ = idx

                if not pd.isnull(intvs_b[i]):
                    ovr_a = merge.index[intvs_a.overlaps(intvs_b[i])]
                    weight_a, sum_a = weighted_dist(ovr_a,idx, False)
                else:
                    weight_a = sum_a = 0

                if not pd.isnull(intvs_a[i]):
                    ovr_b = merge.index[intvs_b.overlaps(intvs_a[i])]
                    weight_b, sum_b = weighted_dist(ovr_b,idx, True)
                else:
                    weight_b = sum_b = 0

                if weight_a + weight_b > 0:
                    mean_ref_dist[idx] = (sum_a+sum_b) / (weight_a+weight_b)

            df.loc[merge.index, "mean_ref_dist"] = mean_ref_dist


    @property
    def read_ids(self):
        return self.alignments["read_id"].unique()

    @property
    def layer_aln_ids(self):
        return self.layers.index.get_level_values("aln_id")

    @property
    def layer_mrefs(self):
        mrefs = self.coords.ref_to_mref(self.layer_refs)
        if self.coords.stranded:
            return mrefs
        else:
            print("UNSTRAND")
            return mrefs[0].union(mrefs[1])
        #return self.layers.index.get_level_values("mref")

    @property
    def layer_refs(self):
        #return self.coords.mref_to_ref_index(self.layer_mrefs)
        return self.layers.index.get_level_values("ref")
