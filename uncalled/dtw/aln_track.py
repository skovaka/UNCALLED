#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
import re
import time
from typing import NamedTuple
import pandas as pd
import copy

from _uncalled import Compare

from ..pore_model import PoreModel
from ..pafstats import parse_paf, PafEntry
from ..argparse import Opt, ref_coords
from .. import config, index 
from ..index import load_index, RefCoord
from .moves import Bcaln

from .layers import LAYER_META


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

    def _init_new(self, track_id, name, desc, conf, model=None, fast5s=None):
        self.id = track_id
        self.name = name
        self.desc = desc
        self.conf = conf
        self.coords = None

        self.fast5s = fast5s #TODO get rid of this

        if model is not None:
            self.model = model 
        elif len(conf.pore_model.name) > 0:
            self.model = PoreModel(conf.pore_model)
        else:
            self.model = None

    def _init_slice(self, p, coords=None, alignments=None, layers=None, order=["fwd", "ref_start"]):
        self._init_new(p.id, p.name, p.desc, p.conf, p.model, p.fast5s)
        self.set_data(coords, alignments, layers, order)

    def _parse_layers(self, df):
        if df.index.names[0] == "pac":
            df = df.rename(index=self.coords.pac_to_ref, level=0)
            df.index.names = ("ref", "aln_id")

            refs = df.index.get_level_values(0)
            if len(df) > 1 and refs[0] > refs[1]:
                df = df.iloc[::-1]

        if self.conf.tracks.mask_indels is not None and ("moves","indel") in df.columns:
            df = df[df["moves","indel"].abs() < self.conf.tracks.mask_indels]


        if self.conf.tracks.mask_skips is not None and ("dtw","events") in df.columns:
            skips = df["dtw","events"] < 1
            if self.conf.tracks.mask_skips == "all":
                df = df[~skips]

            elif self.conf.tracks.mask_skips == "keep_best":
                sdf = df[skips].set_index(("dtw","start"), append=True)
                #sdf.index.names = ("ref", "aln_id", "start")
                sdf["diff"] = (sdf["dtw","current"] - self.model[sdf["dtw","kmer"]]).abs()
                grp = sdf.groupby(level=2, group_keys=False)
                keep_idx = grp["diff"].idxmin()
                sdf = sdf.drop(keep_idx).reset_index(level=2)
                df = df.drop(sdf.index)
            else:
                raise ValueError("Unknown masking mode: " + self.conf.tracks.mask_skips)
            
        return df

    def set_data(self, coords, alignments, layers, order=["fwd", "ref_start"]):

        self.coords = coords
        self.alignments = alignments
        self.layers = self._parse_layers(layers)

        if not self.coords.stranded and (self.all_fwd or self.all_rev):
            self.coords = self.coords.ref_slice(fwd=self.all_fwd)

        isnone = [coords is None, alignments is None, layers is None]
        if np.all(isnone) or len(alignments) == 0:
            return
        elif np.any(isnone):
            raise ValueError("Must specify AlnTrack coords, alignments, and layers")
        self.alignments = self.alignments.sort_values(order)


        #self.layers = self.layers.sort_index()

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


        #layers = self.layers.loc[layer_refs.isin(coords.refs)]
        layers = self.layers.loc[(layer_refs >= coords.refs.min()) & (layer_refs <= coords.refs.max())]

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

        layers["aln_id"] = aln_id
        df = layers.set_index("aln_id", append=True).reorder_levels(["ref","aln_id"])

        meta = LAYER_META.loc[group]
        df = df[meta.index.intersection(df.columns)]
        df = df.astype(meta.loc[df.columns, "dtype"], copy=False)

        df = self._parse_layers(pd.concat({group : df}, names=["group", "layer"], axis=1))
        pd.set_option('display.max_rows', 50)


        if self.layers is None or overwrite:
            self.layers = df #pd.DataFrame({
            #    layer : df[layer].astype(LAYER_META.loc[layer,"dtype"]) 
            #    for layer in df}, columns=df.columns)
        else:
            #TODO don't always parse every layer twice
            #self.layers = pd.concat([self.layers, self._parse_layers(df)], axis=1)
            self.layers = pd.concat([self.layers, df], axis=1)
            #for layer in df:
            #    self.layers[layer] = df[layer]#.astype(LAYER_META.loc[layer,"dtype"])

        return df

    def aln_ref_coord(self, aln_id):
        return RefCoord(*self.alignments[["ref_name","ref_start","ref_end","fwd"]].loc[aln_id])

    def has_group(self, group):
        return group in self.layers.columns.get_level_values(0)


    def calc_layers(self, layers):
        for group, layer in layers:
            if group in self.layers and (group, layer) not in self.layers.columns:
                meta = LAYER_META.loc[(group,layer)]

                #Make sure layer dependencies exist
                if not self.empty and (meta["deps"] is None or len(self.layers.columns.intersection(meta["deps"])) == len(meta["deps"])):

                    fn = meta["fn"]
                    if fn is None:
                        continue
                        #raise ValueError(f"Layer not found: {group}.{layer}")
                    vals = fn(self)
                    self.layers[group,layer] = vals

    def load_mat(self):
        df = self.layers.copy()
        #df["aln_id"] = df.index.get_level_values("aln_id")
        df = df.reset_index()

        self.mat = df.pivot(index="aln_id", columns=["ref"]) \
                     .rename_axis(("group","layer","ref"), axis=1) \
                     .reindex(self.coords.refs, axis=1, level=2) \
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

    def cmp(self, other, calc_jaccard, calc_dist):
        groups_b = other.alignments.groupby("read_id")

        df = pd.DataFrame(
            columns=["aln_b", "group_b", "dist", "jaccard"],
            index = self.layers.index
        )

        for id_a, aln_a in self.alignments.iterrows():
            read_id = aln_a["read_id"]
            self.calc_layers([("dtw","end")])
            dtw_a = self.layers.loc[(slice(None),id_a),"dtw"][["start","end"]]

            for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                self._compare_alns(dtw_a, other, id_b, "dtw", df, calc_jaccard, calc_dist)

        df["group_b"] = "dtw"
        return df#.set_index(["aln_b", "group_b"], append=True)

    def mvcmp(self, other, calc_jaccard, calc_dist):
        if other != self:
            groups_b = other.alignments.groupby("read_id")

        df = pd.DataFrame(
            columns=["aln_a", "aln_b", "group_b", "dist", "jaccard"],
            index = self.layer_refs
        )

        for id_a, aln_a in self.alignments.iterrows():
            self.calc_layers([("dtw","end")])
            dtw = self.layers.loc[(slice(None),id_a),"dtw"][["start","end"]]

            if dtw is None:
                continue

            if other == self:
                self._compare_alns(dtw, self, id_a, "moves", df, calc_jaccard, calc_dist)
            else:
                read_id = aln_a["read_id"]
                for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                    self._compare_alns(dtw, other, id_b, "moves", df)

        df["group_b"] = "moves"
        return df#.set_index(["aln_b", "group_b"], append=True)

    def _compare_alns(self, aln_a, other, id_b, group, df, calc_jaccard=True, calc_dist=True):

        other.calc_layers([(group,"end")])

        aln_b = other.layers \
                   .loc[(slice(None),id_b),group][["start","end"]] \
                   .reset_index(level="aln_id") \
                   .rename(columns={"aln_id" : "aln_b"})
                

        alns_a = aln_a.index.unique(1)
        if len(alns_a) != 1:
            raise ValueError("Can only compare two alignments at a time")
        
        has_pac = "pac" in aln_b.index.name
        flip = self.all_fwd == self.conf.is_rna

        def coords(df, track):
            df = df.dropna().sort_index().reset_index()
            if flip:
                #df = df[::-1]
                #df["ref"] = -df["ref"]
                end = -df["start"]
                df["start"] = -df["end"]
                df["end"] = end
            return AlnCoords(df)


        coords_a = coords(aln_a, self)
        coords_b = coords(aln_b, other)

        compare = Compare(coords_a, coords_b)
        #cmp_df = compare.to_numpy()
        #idx = (cmp_df["ref"],slice(None))
        #df.loc[idx,"jaccard"] = cmp_df["jaccard"]
        #df.loc[idx,"dist"] = cmp_df["dist"]

        cmp_df = pd.DataFrame(compare.to_numpy()).dropna(how="all")

        #cmp_df["aln_b"] = alns_b[0]
        if has_pac:
            cmp_df["pac"] = self.coords.ref_to_pac(pd.Index(cmp_df["ref"]))
            cmp_df = cmp_df.set_index("pac")
        else:
            cmp_df = cmp_df.set_index("ref")
        df["jaccard"] = cmp_df["jaccard"]
        df["dist"] = cmp_df["dist"]
        df["aln_a"] = alns_a[0]
        df["aln_b"] = id_b



    @property
    def read_ids(self):
        return self.alignments["read_id"].unique()

    @property
    def layer_aln_ids(self):
        return self.layers.index.get_level_values("aln_id")

    @property
    def layer_pacs(self):
        return self.coords.ref_to_pac(self.layer_refs)

    @property
    def layer_fwds(self):
        return self.alignments.loc[self.layer_aln_ids, "fwd"]

    @property
    def layer_strands(self):
        return self.layer_fwds.map({True : "+", False : "-"})

    @property
    def layers_pac_index(self):
        index = pd.MultiIndex.from_arrays([self.coords.ref_to_pac(self.layer_refs), self.layers.index.get_level_values("aln_id")], names=["pac","aln_id"])
        return self.layers.set_index(index, drop=True)

    @property
    def layers_desc_index(self):
        index = pd.MultiIndex.from_arrays(
            [self.layer_refs, self.layer_strands, self.layer_aln_ids])
        return pd.concat({self.coords.ref_name : self.layers.set_index(index, drop=True)}, names=[("ref","name"),("ref","coord"),("ref","strand"),"aln_id"])

        stats.index = pd.MultiIndex.from_product([
            [chunk.coords.ref_name], stats.index, ["+" if chunk.coords.fwd else "-"]
        ])

    @property
    def layer_refs(self):
        return self.layers.index.get_level_values("ref")
