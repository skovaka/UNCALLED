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

from .read_aln import ReadAln, RefCoord, LayerMeta, LAYER_META
from .track_io import TrackSQL

from ..pafstats import parse_paf, PafEntry
from ..argparse import Opt, ref_coords
from .. import nt, PoreModel, config, index 
from ..index import load_index

KMER_LAYER       = "kmer"
CURRENT_LAYER    = "current"
DWELL_LAYER      = "dwell"
MODEL_DIFF_LAYER = "model_diff"
DEFAULT_LAYERS = [CURRENT_LAYER, DWELL_LAYER, MODEL_DIFF_LAYER]

LAYER_META.update({
    "dwell" : LayerMeta(float, "Dwell Time (ms/nt)"),
    "model_diff" : LayerMeta(float, "Model pA Difference"),
    "bcerr"      : LayerMeta(float, "BC Alignment Errors")
})

class AlnTrackParams(config.ParamGroup):
    _name = "track"
AlnTrackParams._def_params(
    ("path", None, None, "Path to directory where alignments are stored"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("name", None, str, "Short unique identifier for the track"),
    ("index_prefix", None, str, "BWA index prefix"),
    ("load_mat", True, bool, "If true will load a matrix containing specified layers from all reads overlapping reference bounds"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("layers", DEFAULT_LAYERS, list, "Layers to load"),
    ("mode", "r", str, "Read (r) or write (w) mode"),
    ("overwrite", False, bool, "Will overwrite existing directories if True"),
    ignore_toml={"mode", "overwrite", "load_mat"}
)


class AlnTrack:
    CONF_FNAME = "conf.toml"
    ALN_DIR = "alns"

    #TODO make fast5_file alias for filename in Fast5Reader
    #also, in fast5 reader make static fast5 filename parser

    WRITE_MODE = "w"
    READ_MODE = "r"
    MODES = {WRITE_MODE, READ_MODE}

    CMP_LAYERS = [CURRENT_LAYER, DWELL_LAYER]

    def get_bcerr_layer(self, aln):
        bcerr = aln.bcerr#.reindex(aln.aln.index)
        ret = pd.Series(np.nan, bcerr.index)
        subs = bcerr[bcerr["type"]=="SUB"]
        ret[subs.index] = subs["seq"].replace({"A":0,"C":1,"G":2,"T":3})
        ret[bcerr["type"]=="DEL"] = 4
        ret[bcerr["type"]=="INS"] = 5
        return ret

    LAYER_FNS = {
        #"id" : (
        #    lambda self,a: a.id),
        "kmer" : (
            lambda self,a: self.load_aln_kmers(a)),
        "current" : (
            lambda self,a: a.aln["current"]),
        "dwell" : (
            lambda self,a: 1000 * a.aln['length'] / self.conf.read_buffer.sample_rate),
        "model_diff" : (
            lambda self,a: a.aln["current"] - self.model[a.aln["kmer"]]),
        "bcerr" : get_bcerr_layer,
    }

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("track", *args, **kwargs)

        self.in_mem = self.prms.path is None

        self.read_ids = set()

        if self.prms.name is None:
            self.prms.name = os.path.basename(self.prms.path)

        if not self.in_mem:
            self.db = TrackSQL(self.prms.path)

            if self.prms.mode == self.READ_MODE:
                self.id, self.desc, toml, groups = self.db.query_track(self.prms.name)
                self.conf.load_toml(text=toml)

        self.index = load_index(self.prms.index_prefix)
        self.model = PoreModel(self.conf.pore_model)
        self.read_aln = None
        self.mat = None
        self.coords = None

        self.set_ref_bounds(self.prms.ref_bounds)

        if self.conf.dtw.mm2_paf is not None:
            read_filter = set(self.conf.fast5_reader.read_filter)

        if not self.in_mem:
            if self.prms.load_mat and self.prms.ref_bounds is not None:
                self.load_region(self.prms.ref_bounds)

            elif self.prms.mode == self.WRITE_MODE:
                self.prev_fast5 = (None, None)
                self.prev_read = None
                self.prev_aln = -1

                self.db.init_tables()
                self.id = self.db.init_track(self)

    def __contains__(self, read_id):
        return read_id in self.read_ids

    def _load_fast5_index(self):
        self.read_ids = set(self.fname_mapping.index)
        if len(self.conf.fast5_reader.read_filter) > 0:
            self.read_ids = self.read_ids & set(self.conf.fast5_reader.read_filter)

    #TODO move to RefIndex python wrapper? 
    def _ref_coords_to_mrefs(self, ref_coords, fwd=None):
        if fwd is None:
            fwd = ref_coords.fwd
        start, end = self.index.ref_to_mref(
            ref_coords.name, ref_coords.start, ref_coords.end-nt.K+1, fwd, self.conf.is_rna)
        return pd.RangeIndex(start, end).rename("mref")


    def set_ref_bounds(self, ref_bounds):

        if ref_bounds == None:
            return

        self.coords = self.index.get_coord_space(ref_bounds, self.conf.is_rna)

        self.width = len(self.coords)
        self.height = None

    def _group_layers(self, name, df):
        return pd.concat({name : df}, names=["group", "layer"], axis=1)
        

    def init_alignment(self, read, group_name, layers):

        #TODO do this here instead of in bcaln or whatever
        #if self.coords is not None:
        #    mrefs = self.coords.mrefs.intersection(layers.index)
        #    if len(mrefs) == 0: 
        #        return None
        #    layers = layers.reindex(index=mrefs)

        fast5 = read.f5.filename
        read_id = read.id

        print(layers)
        mref_start = layers.index.min()
        mref_end = layers.index.max()
        samp_start = layers["sample"].min()
        samp_end = layers["sample"].max()

        ref_bounds = self.index.mref_to_ref_bound(mref_start, mref_end, not self.conf.is_rna)

        self.prev_aln += 1
        aln_id = self.prev_aln

        self.alignments = pd.DataFrame({
                "id" : [aln_id],
                "track_id" : [self.id],
                "read_id" : [read_id],
                "ref_name" : [ref_bounds.ref_name],
                "ref_start" : [ref_bounds.start],
                "ref_end" : [ref_bounds.end],
                "fwd" :     [ref_bounds.fwd],
                "samp_start" : [samp_start],
                "samp_end" : [samp_end],
                "tags" : [""]}).set_index("id")

        self.layers = self._group_layers(group_name, layers)
        self.layers.index = pd.MultiIndex.from_product([self.layers.index, [aln_id]], names=["mref", "aln_id"])
        print(self.layers)

        #TODO iterate hierarchically through fast5s/reads
        if fast5 == self.prev_fast5[0]:
            fast5_id = self.prev_fast5[1]
        else:
            fast5_id = self.db.init_fast5(fast5)
            self.prev_fast5 = (fast5, fast5_id)

        if self.prev_read != read_id:
            self.db.init_read(read_id, fast5_id)
            self.prev_read = read_id

        self.db.init_alignment(self.alignments)
        self.db.write_layers(group_name, self.layers[group_name])

        return aln_id

    def load_aln_kmers(self, aln=None, store=True):
        if aln is None:
            aln = self.read_aln

        kmers = pd.Series(
            self.index.get_kmers(aln.mref_start-nt.K+1, aln.mref_end, aln.is_rna),
            #aln.mrefs
            pd.RangeIndex(aln.mref_start, aln.mref_end)
        )

        if store:
            self.read_aln.dtw["kmer"] = kmers
        return kmers

    def save_aln(self, aln=None):
        if self.prms.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")
        
        if aln is not None:
            self.read_aln = aln

        self.db.write_layers("dtw", self.read_aln.dtw)


    def load_read(self, read_id=None, load_kmers=True):
        if self.read_aln is not None and read_id == self.read_aln.read_id: 
            return self.read_aln

        aln, dtw = self.db.query_read(read_id, self.coords)

        self.read_aln = ReadAln(self.id, read_id, index=self.index, is_rna=self.conf.is_rna)
        self.read_aln.set_dtw(dtw)

        if not self.read_aln.empty:
            if load_kmers:
                self.load_aln_kmers()

            for layer in self.prms.layers:
                if not layer in self.read_aln.dtw.columns:
                    self.read_aln.dtw[layer] = self.LAYER_FNS[layer](self, self.read_aln)

        return self.read_aln

    def set_layers(self, layers):
        if layers is not None:
            self.prms.layers = layers
        self.layer_idxs = {layer : i for i,layer in enumerate(self.prms.layers)}

    #TODO parse mm2 every time to enable changing bounds
    #eventually use some kind of tabix-like indexing
    def load_region(self, ref_bounds=None, layers=None):
        self.set_layers(layers)
        if ref_bounds is not None:
            self.prms.ref_bounds = ref_bounds

        read_meta = defaultdict(list)

        read_rows = defaultdict(list)
        mask_rows = list()

        ##TODO make index take RefCoord (combine with RefLoc)

        self.set_ref_bounds(self.prms.ref_bounds)

        self.alignments, self.df = self.db.query_alns(self.id, self.coords, self.prms.full_overlap)

        self.df.set_index(["mref", "aln_id"], inplace=True)
        self.alignments.sort_values("ref_start", inplace=True)

        mat_index = pd.MultiIndex.from_product([["current", "start", "length"], self.coords.refs])

        self.mat = self.df.reset_index().pivot(index="aln_id", columns="mref") \
                   .rename(columns=self.coords.mref_to_ref, level=1) \
                   .rename_axis(("layer","ref"), axis=1) \
                   .reindex(mat_index, axis=1)

        self.mat = self.mat.reindex(self.alignments["id"], copy=False)

        self.has_fwd = np.any(self.alignments['fwd'])
        self.has_rev = not np.all(self.alignments['fwd'])

        self.height = len(self.alignments)


        return self.mat

    def add_layer(self, name, mat):
        self.mat = np.ma.concatenate([self.mat, [mat]])
        self.set_layers(self.prms.layers + [name])

    def get_pileup(self, layer):
        return np.flip(np.sort(self[layer], axis=0), axis=0)

    def close(self):
        if self.in_mem: return
        self.db.close()

    #@property
    #def name(self):
    #    if self.prms.path is None:
    #        return self.fast5s
    #    return self.prms.path.split("/")[-1]

    @property
    def read_count(self):
        return len(self.fname_mapping)
    
    def sort_coord(self, layer, ref):
        #if isinstance(layer, str):
        #    layer = self.layer_idxs[layer]
        order = (-self.mat[layer,ref].fillna(0)).argsort()
        self.sort(order)

    def sort(self, order):
        self.mat = self.mat.iloc[order]
        self.alignments = self.alignments.iloc[order]

    @property
    def ref_id(self):
        return self.index.get_ref_id(self.ref_name)

    @property
    def ref_name(self):
        return self.prms.ref_bounds.name

    @property
    def ref_start(self):
        return self.prms.ref_bounds.start

    @property
    def ref_end(self):
        return self.prms.ref_bounds.end

    def calc_ks(self, track_b):
        ks_stats = np.zeros((len(self.CMP_LAYERS), self.width))

        for i,l in enumerate(AlnTrack.CMP_LAYERS):
            for j,rf in enumerate(self.coords.mrefs[True]):
                a = self[l,:,rf]
                b = track_b[l,:,rf]
                ks = scipy.stats.mstats.ks_2samp(a,b,mode="asymp")
                ks_stats[i][j] = ks[0]

        return ks_stats

    def calc_pca(self, layer, ref_start, ref_end, n_components=2):
        x = self[layer,:,ref_start:ref_end].T
        pc = PCA(n_components=n_components).fit_transform(x)
        data = {"read_id" : self.alignments["id"]}
        for c in range(n_components):
            data["pc%d" % (c+1)] = pc[:,c]
        return pd.DataFrame(data).set_index("read_id")
    
    
    def sort_pca(self, layer, ref_start, ref_end):
        pc = self.calc_pca(layer, ref_start, ref_end, 1)
        self.sort(pc["pc1"].argsort())


    def __getitem__(self, key):
        if isinstance(key, str):
            #key = self.layer_idxs[key]
            return self.mat[key]
        if len(key) == 3:
            return self.mat[key[0],key[2]][key[1]]
        return self.mat[key[0]][key[1]]
    
    @property
    def kmer(self):
        return self.mat[self.KMER_LAYER]

    @property
    def pa(self):
        return self.mat[self.PA_LAYER]
    
    @property
    def dwell(self):
        return self.mat[self.DWELL_LAYER]
    
    @property
    def pa_diff(self):
        return self.mat[self.PA_DIFF_LAYER]

#TODO turn into "AlnTrackIO"
def _load_tracks(tracks, conf, force_list=False):

    is_list = isinstance(tracks, list)
    if not is_list:
        tracks = [tracks]
    else:
        ret = list()
    
    for track in tracks:
        if isinstance(track, str):
            track = AlnTrack(track, conf=conf)
        elif not isinstance(track, AlnTrack):
            raise RuntimeError("AlnTrack must either be path to alignment track or AlnTrack instance")
        
        if is_list or force_list:
            ret.append(track)
        else:
            return track
    return ret

def _load_track_arg(arg_track, conf_track, conf):                                       
    track = arg_track if arg_track is not None else conf_track                          
    if isinstance(track, str):
        return AlnTrack(track, conf=conf)                                                  
    elif isinstance(track, AlnTrack):
        return track
    raise RuntimeError("AlnTrack must either be path to alignment track or AlnTrack instance")

class CompareParams(config.ParamGroup):
    _name = "compare"
CompareParams._def_params(
    ("track_a", None, None, "DTW AlnTrack A"),
    ("track_b", None, None, "DTW AlnTrack B"),
    ("subcmd", "ks", str, "Analysis to perform")
)

COMPARE_OPTS = (
    Opt("ref_bounds", "track", type=ref_coords),
    Opt("track_a", type=str),
    Opt("track_b", type=str, nargs="?"),
    Opt(("-f", "--full-overlap"), "track", action="store_true"),
    #Opt(("-o", "--outfile"), type=str, default=None),
)

#def compare(track_a=None, track_b=None, full_overlap=None, conf=None):
def compare(*args, **kwargs):
    """Outputs a TSV file conaining Kolmogorov–Smirnov test statistics comparing the current and dwell time of two alignment tracks"""
    conf, prms = config._init_group("compare", *args, **kwargs)

    track_a = _load_track_arg(None, prms.track_a, conf)
    track_b = _load_track_arg(None, prms.track_b, conf)

    ks_stats = np.recarray(
        track_a.width,
        dtype=[(l, '<f8') for l in AlnTrack.CMP_LAYERS]
    )

    for i,l in enumerate(AlnTrack.CMP_LAYERS):
        for j,rf in enumerate(track_a.ref_coords.index):
            a = track_a[l,:,rf]
            b = track_b[l,:,rf]
            ks = scipy.stats.mstats.ks_2samp(a,b,mode="asymp")
            ks_stats[j][i] = ks[0]

    return pd.DataFrame(
        ks_stats, 
        index = track_a.ref_coords.index.rename("ref")
    )


METHOD_COMPARE_OPTS = (
    Opt("track_a", "browser"),
    Opt("track_b", "browser"),
    #Opt(("-R", "--ref-bounds"), "dtw", type=ref_coords),
    Opt(("-f", "--full-overlap"), "browser", action="store_true"),
)

def method_compare(track_a=None, track_b=None, full_overlap=None, conf=None):
    if conf is None:
        conf = conf if conf is not None else config.Config()
    if full_overlap is not None:
        conf.browser.full_overlap = full_overlap
    track_a = _load_track_arg(track_a, conf.browser.track_a, conf)
    track_b = _load_track_arg(track_b, conf.browser.track_b, conf)

    common_reads = set(track_a.read_ids) & set(track_b.read_ids)
    if len(common_reads) == 0:
        sys.stderr.write("Error: method_compare tracks must have reads in common\n")
        return

    for read in common_reads:
        aln_a = track_a.load_read(read)
        aln_b = track_b.load_read(read)

def method_compare_aln(aln_a, aln_b):
    merge = aln_a.aln.join(aln_b.aln, lsuffix="_a", rsuffix="_b").dropna().set_index("mref_a")

    def get_ends(suff):
        return merge["start_"+suff]+merge["length_"+suff]

    st_a = merge["start_a"]
    st_b = merge["start_b"]
    en_a = get_ends("a")
    en_b = get_ends("b")

    mid_a = st_a + merge["length_a"] 
    mid_b = st_b + merge["length_b"] 

    st_min = np.minimum(st_a, st_b)
    st_max = np.maximum(st_a, st_b)
    en_min = np.minimum(en_a, en_b)
    en_max = np.maximum(en_a, en_b)

    df = pd.DataFrame({
        "jaccard" : np.maximum(0, en_min-st_max)/(en_max-st_min),
        "centroid_diff" : mid_a-mid_b,
        "dwell_diff" : merge["length_a"]-merge["length_b"],
    })
    return df
        
#class AlnTrackParams(config.ParamGroup):
#    _name = "track"
#BrowserParams._def_params(
#    ("track_a", None, str, "Path to directory where alignments are stored"),
#    ("track_b", None, str, "Path to directory where alignments are stored"),
#    ("stats", None, str, "Path to directory where alignments are stored"),

REFSTATS_OPTS = (
    Opt("ref_bounds", "track", type=ref_coords),
    Opt("track_in", type=str),
    Opt(("-f", "--full-overlap"), "track", action="store_true"),
    Opt(("-m", "--pore-model"), "mapper", default=None),
    Opt("--rna", fn="set_r94_rna"),
)

def refstats(
        track=None, 
        full_overlap=None, 
        layers=[CURRENT_LAYER, DWELL_LAYER, MODEL_DIFF_LAYER], #TODO use strings
        stats=["mean", "variance", "skewness", "kurtosis"], 
        conf=None):

    if conf is None:
        conf = conf if conf is not None else config.Config()
    if full_overlap is not None:
        conf.browser.full_overlap = full_overlap
    track = _load_track_arg(track, getattr(conf, "track_in", None), conf)

    vals = dict()
    data = list()
    names = list()

    mask = ~track[0].mask
    vals["cov"] = np.sum(mask, axis=0)
    vals["kmer"] = [
        nt.kmer_to_str(int(k)) 
        for k in np.ma.median(track[KMER_LAYER], axis=0).data
    ]

    for l in layers:
        desc = scipy.stats.mstats.describe(track[l], axis=0)
        for stat in stats:
            vals["%s_%s" % (l, stat)] = getattr(desc, stat)

    return pd.DataFrame(vals, index=track.ref_coords.index.rename("ref"))
