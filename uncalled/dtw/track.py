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
import scipy.stats
import copy

import matplotlib.pyplot as plt

from . import ReadAln, RefCoord

from ..pafstats import parse_paf, PafEntry
from .. import config 
from ..config import Config, Opt
from .. import BwaIndex, nt, PoreModel

KMER_LAYER       = "kmer"
CURRENT_LAYER    = "current"
DWELL_LAYER      = "dwell"
MODEL_DIFF_LAYER = "model_diff"
DEFAULT_LAYERS = [KMER_LAYER, CURRENT_LAYER, DWELL_LAYER, MODEL_DIFF_LAYER]
 
def ref_coords(coord_str):             
    spl = coord_str.split(":")         
    ch = spl[0]                        
    st,en = spl[1].split("-")          
                                       
    coord = (ch, int(st), int(en))     
                                       
    if len(spl) == 2:                  
        return coord                   
    else:                              
        return coord + (spl[2] == "+",)

class TrackParams(config.ParamGroup):
    _name = "track"
TrackParams._def_params(
    ("path", None, str, "Path to directory where alignments are stored"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("load_mat", True, bool, "If true will load a matrix containing specified layers from all reads overlapping reference bounds"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("layers", DEFAULT_LAYERS, list, "Layers to load"),
    ("mode", "r", str, "Read (r) or write (w) mode"),
    ("overwrite", False, bool, "Will overwrite existing directories if True"),
    ignore_toml={"mode", "overwrite"}
)

class Track:
    CONF_FNAME = "conf.toml"
    ALN_DIR = "alns"
    ALN_SUFFIX = ".pkl"

    INDEX_FNAME = "filename_mapping.txt"
    INDEX_HEADER = "read_id\tfilename\taln_file"

    #TODO make fast5_file alias for filename in Fast5Reader
    #also, in fast5 reader make static fast5 filename parser

    WRITE_MODE = "w"
    READ_MODE = "r"
    MODES = {WRITE_MODE, READ_MODE}

    CMP_LAYERS = [CURRENT_LAYER, DWELL_LAYER]

    LAYER_META = [
        ("K-mer",              False),
        ("Current (pA)",       True),
        ("Dwell Time (ms/bp)", True),
        ("pA Difference",      True)
    ]

    LAYER_FNS = {
        KMER_LAYER : (
            lambda self,df: df["kmer"]),
        CURRENT_LAYER : (
            lambda self,df: df["mean"]),
        DWELL_LAYER : (
            lambda self,df: df['length'] / self.conf.read_buffer.sample_rate),
        MODEL_DIFF_LAYER : (
            lambda self,df: df["mean"] - self.model[df["kmer"]]),
    }

    #def __init__(self, path, mode="r", ref_bounds=None, full_overlap=None, conf=None, overwrite=False, index=None, mm2s=None):
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("track", *args, **kwargs)

        if self.prms.mode == self.WRITE_MODE:
            os.makedirs(self.aln_dir, exist_ok=self.prms.overwrite)

        self.fname_mapping_file = open(self.fname_mapping_filename, self.prms.mode)

        if self.prms.mode == self.READ_MODE:
            self.conf.load_toml(self.config_filename)

            if len(self.conf.fast5_reader.fast5_index) == 0:
                self.conf.fast5_reader.fast5_index = self.fname_mapping_filename

            self._load_index()

        elif self.prms.mode == self.WRITE_MODE:
            self.conf.to_toml(self.config_filename)
            self.fname_mapping_file.write(self.INDEX_HEADER + "\n")

        #TODO arguments overload conf params
        if self.conf.align.mm2_paf is not None:
            read_filter = set(self.conf.fast5_reader.read_filter)
            self.mm2s = {p.qr_name : p
                     for p in parse_paf(
                        self.conf.align.mm2_paf,
                        ref_bounds=self.prms.ref_bounds,
                        full_overlap=self.prms.full_overlap,
            )}

        #TODO static bwa_index parameters, or instance
        if len(self.conf.mapper.bwa_prefix) > 0:
            self.index = BwaIndex(self.conf.mapper.bwa_prefix, True)

        self.model = PoreModel(self.conf.pore_model)

        if self.prms.load_mat and self.prms.ref_bounds is not None:
            self.load_region(self.prms.ref_bounds)

    #@property
    #def read_ids(self):
    #    return list(self.fname_mapping.index)

    def __contains__(self, read_id):
        return read_id in self.read_ids

    def _load_index(self):
        self.fname_mapping = pd.read_csv(self.fname_mapping_file, sep="\t", index_col="read_id")
        self.read_ids = set(self.fname_mapping.index)

    def save_aln(self, aln, fast5_fname):
        if self.prms.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")

        aln_fname = self.aln_fname(aln.read_id)
        aln.df.sort_index().to_pickle(aln_fname)
        self.fname_mapping_file.write("\t".join([aln.read_id, fast5_fname, aln_fname]) + "\n")

    def load_aln(self, read_id, ref_bounds=None):
        mm2 = self.mm2s[read_id]
        df = pd.read_pickle(self.aln_fname(read_id)).sort_index()
        return ReadAln(self.index, mm2, df, is_rna=not self.conf.read_buffer.seq_fwd, ref_bounds=ref_bounds)

    def set_layers(self, layers):
        if layers is not None:
            self.prms.layers = layers
        self.layer_idxs = {layer : i for i,layer in enumerate(self.prms.layers)}

    #TODO parse mm2 every time to enable changing bounds
    #eventually use some kind of tabix-like indexing
    def load_region(self, ref_bounds=None, layers=None):

        #self.mat = TrackMatrix(self, ref_bounds, self.mm2s, conf=self.conf)

        self.set_layers(layers)
        if ref_bounds is not None:
            self.prms.ref_bounds = ref_bounds

        self.width = self.ref_end-self.ref_start
        self.height = None

        read_meta = defaultdict(list)

        read_rows = defaultdict(list)
        mask_rows = list()

        self.ref_coords = pd.Series(
                np.arange(self.width, dtype=int),
                index=pd.RangeIndex(self.ref_start, self.ref_end)
        )

        def _add_row(df, mm2_paf):
            dtw_roi = df.loc[self.ref_start:self.ref_end-1]

            xs = self.ref_coords.reindex(
                self.ref_coords.index.intersection(dtw_roi.index)
            )

            roi_mask = np.ones(self.width)#, dtype=bool)
            roi_mask[xs] = False
            mask_rows.append(roi_mask)

            #pa_diffs = dtw_roi["mean"] - self.model[dtw_roi["kmer"]]
            #pa_diffs = self.model.abs_diff(dtw_roi['mean'], dtw_roi['kmer'])
            #dwell = 1000 * dtw_roi['length'] / self.conf.read_buffer.sample_rate

            #def _add_layer_row(layer, vals):
            for layer in self.prms.layers:
                row = np.zeros(self.width)
                row[xs] = self.LAYER_FNS[layer](self, dtw_roi)
                read_rows[layer].append(row)

            #_add_layer_row(self.KMER_LAYER, dtw_roi['kmer'])
            #_add_layer_row(self.PA_LAYER, dtw_roi['mean'])
            #_add_layer_row(self.PA_DIFF_LAYER, pa_diffs)
            #_add_layer_row(self.DWELL_LAYER, dwell)

            read_meta['ref_start'].append(df.index.min())
            read_meta['id'].append(mm2_paf.qr_name)
            read_meta['fwd'].append(mm2_paf.is_fwd)

            #self.mm2s[mm2_paf.qr_name] = mm2_paf

        for read_id,read in self.fname_mapping.iterrows():
            mm2 = self.mm2s.get(read_id, None)
            if mm2 is None: 
                continue
            df = pd.read_pickle(read["aln_file"])
            df.sort_index(inplace=True)
            _add_row(df, mm2)

        self.reads = pd.DataFrame(read_meta) \
                     .sort_values(['fwd', 'ref_start'], ascending=[False, True])

        self.has_fwd = np.any(self.reads['fwd'])
        self.has_rev = not np.all(self.reads['fwd'])

        read_order = self.reads.index.to_numpy()

        mat = np.stack([
            np.stack(read_rows[l])[read_order] for l in self.prms.layers
        ])
        mask = np.stack([np.stack(mask_rows)[read_order].astype(bool) for _ in self.prms.layers])

        self.mat = np.ma.masked_array(mat, mask=mask)

        self.height = len(self.reads)

        #self.layer_extrema = np.array([
        #    [np.min(layer), np.max(layer)]
        #    for layer in self.mat
        #])
        #TODO define get_soft_minmax(num_stdvs):

        self.norms = [Normalize(np.min(layer), np.max(layer)) for layer in self.mat]
        for l in [CURRENT_LAYER, DWELL_LAYER]:
            i = self.layer_idxs[l]
            layer = self.mat[i]
            self.norms[i].vmax = min(
                layer.max(),
                np.ma.median(layer) + 2 * layer.std()
            )

        return mat


    def close(self):
        self.fname_mapping_file.close()

    @property
    def config_filename(self):
        return os.path.join(self.prms.path, self.CONF_FNAME)

    @property
    def fname_mapping_filename(self):
        return os.path.join(self.prms.path, self.INDEX_FNAME)
    
    @property
    def aln_dir(self):
        return os.path.join(self.prms.path, self.ALN_DIR)
    
    @property
    def read_count(self):
        return len(self.fname_mapping)
    
    def aln_fname(self, read_id):
        return os.path.join(self.aln_dir, read_id+self.ALN_SUFFIX)

    def sort(self, layer, ref):
        order = np.argsort(self.mat[layer,:,ref-self.ref_start])
        self.mat = self.mat[:,order,:]
        self.reads = self.reads.iloc[order]

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

        for i,l in enumerate(self.CMP_LAYERS):
            for rf in range(self.width):
                a = self[l,:,rf]
                b = track_b[l,:,rf]
                ks_stats[i,rf] = scipy.stats.ks_2samp(a,b,mode="asymp")[0]

        return ks_stats

    def __getitem__(self, key):
        if isinstance(key, str):
            key = self.layer_idxs[key]
        elif isinstance(key, tuple) and isinstance(key[0], str):
            key = (self.layer_idxs[key[0]],) + key[1:]
        elif isinstance(key, list):
            key = [self.layer_idxs[l] if isinstance(l, str) else l for l in key]
        return self.mat.__getitem__(key)
    
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

def _load_track_arg(arg_track, conf_track, conf):                                       
    track = arg_track if arg_track is not None else conf_track                          
    if isinstance(track, str):                                                          
        return Track(track, conf=conf)                                                  
    elif isinstance(track, Track):                                                      
        return track                                                                    
    raise RuntimeError("Track must either be path to alignment track or Track instance")


COMPARE_OPTS = (
    Opt("ref_bounds", "track", type=ref_coords),
    Opt("track_a", type=str),
    Opt("track_b", type=str, nargs="?"),
    Opt(("-f", "--full-overlap"), "track", action="store_true"),
    #Opt(("-o", "--outfile"), type=str, default=None),
)

def compare(track_a=None, track_b=None, full_overlap=None, conf=None):
    """Outputs a TSV file conaining Kolmogorov–Smirnov test statistics comparing the current and dwell time of two alignment tracks"""

    if conf is None:
        conf = conf if conf is not None else Config()

    if full_overlap is not None:
        conf.browser.full_overlap = full_overlap

    track_a = _load_track_arg(track_a, conf.track_a, conf)
    track_b = _load_track_arg(track_b, conf.track_b, conf)

    ks_stats = np.recarray(
        track_a.width,
        dtype=[(l, '<f8') for l in Track.CMP_LAYERS]
    )

    for i,l in enumerate(Track.CMP_LAYERS):
        for rf in range(track_a.width):
            a = track_a[l,:,rf]
            b = track_b[l,:,rf]
            ks = scipy.stats.mstats.ks_2samp(a,b,mode="asymp")
            ks_stats[rf][i] = ks[0]

    return pd.DataFrame(
        ks_stats, 
        index = track_a.ref_coords.index.rename("ref")
    )


METHOD_COMPARE_OPTS = (
    Opt("track_a", "browser"),
    Opt("track_b", "browser"),
    #Opt(("-R", "--ref-bounds"), "align", type=ref_coords),
    Opt(("-f", "--full-overlap"), "browser", action="store_true"),
)

def method_compare(track_a=None, track_b=None, full_overlap=None, conf=None):
    if conf is None:
        conf = conf if conf is not None else Config()
    if full_overlap is not None:
        conf.browser.full_overlap = full_overlap
    track_a = _load_track_arg(track_a, conf.browser.track_a, conf)
    track_b = _load_track_arg(track_b, conf.browser.track_b, conf)

    common_reads = set(track_a.read_ids) & set(track_b.read_ids)
    if len(common_reads) == 0:
        sys.stderr.write("Error: method_compare tracks must have reads in common\n")
        return

    for read in common_reads:
        aln_a = track_a.load_aln(read)
        aln_b = track_b.load_aln(read)
        print (method_compmare_aln(aln_a, aln_b))

def method_compare_aln(aln_a, aln_b):
    merge = aln_a.df.join(aln_b.df, lsuffix="_a", rsuffix="_b").dropna().set_index("refmir_a")

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
    #print(merge[["start_a","length_a","start_b","length_b"]])
    return df
        
#class TrackParams(config.ParamGroup):
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
        conf = conf if conf is not None else Config()
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
