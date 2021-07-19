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

from . import ReadAln, RefCoord, LayerMeta, LAYER_META

from ..pafstats import parse_paf, PafEntry
from ..config import Config, Opt
from .. import nt, PoreModel, config, index 
from ..index import load_index

KMER_LAYER       = "kmer"
CURRENT_LAYER    = "current"
DWELL_LAYER      = "dwell"
MODEL_DIFF_LAYER = "model_diff"
DEFAULT_LAYERS = [KMER_LAYER, CURRENT_LAYER, DWELL_LAYER, MODEL_DIFF_LAYER]

LAYER_META.update({
    "dwell" : LayerMeta(float, "Dwell Time (ms/nt)"),
    "model_diff" : LayerMeta(float, "Model pA Difference"),
    "bcerr"      : LayerMeta(float, "BC Alignment Errors")
})
 
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
    ("path", None, None, "Path to directory where alignments are stored"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("index_prefix", None, str, "BWA index prefix"),
    ("load_mat", True, bool, "If true will load a matrix containing specified layers from all reads overlapping reference bounds"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("layers", DEFAULT_LAYERS, list, "Layers to load"),
    ("mode", "r", str, "Read (r) or write (w) mode"),
    ("overwrite", False, bool, "Will overwrite existing directories if True"),
    ignore_toml={"mode", "overwrite"}
)

class Track:
    HDF_FNAME = "alns.h5"
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

    #LAYER_FNS = {
    #    "kmer" : (
    #        lambda self,df: df["kmer"]),
    #    "current" : (
    #        lambda self,df: df["current"]),
    #    "dwell" : (
    #        lambda self,df: 1000 * df['length'] / self.conf.read_buffer.sample_rate),
    #    "model_diff" : (
    #        lambda self,df: df["current"] - self.model[df["kmer"]]),
    #}

    def get_bcerr_layer(self, aln):
        bcerr = aln.bcerr.reindex(aln.aln.index)
        ret = pd.Series(np.nan, bcerr.index)
        subs = bcerr[bcerr["type"]=="SUB"]
        ret[subs.index] = subs["seq"].replace({"A":0,"C":1,"G":2,"T":3})
        ret[bcerr["type"]=="DEL"] = 4
        ret[bcerr["type"]=="INS"] = 5
        return ret

    LAYER_FNS = {
        "kmer" : (
            lambda self,a: a.aln["kmer"]),
        "current" : (
            lambda self,a: a.aln["current"]),
        "dwell" : (
            lambda self,a: 1000 * a.aln['length'] / self.conf.read_buffer.sample_rate),
        "model_diff" : (
            lambda self,a: a.aln["current"] - self.model[a.aln["kmer"]]),
        "bcerr" : get_bcerr_layer
    }

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("track", *args, **kwargs)

        self.in_mem = self.prms.path is None

        self.read_ids = set()
        if not self.in_mem:
            if self.prms.mode == self.WRITE_MODE:
                os.makedirs(self.aln_dir, exist_ok=self.prms.overwrite)

            self.fname_mapping_file = open(self.fname_mapping_filename, self.prms.mode)
            self.hdf = pd.HDFStore(self.hdf_filename, mode=self.prms.mode, complib="lzo")

            if self.prms.mode == self.READ_MODE:
                self.conf.load_toml(self.config_filename)

                if len(self.conf.fast5_reader.fast5_index) == 0:
                    self.conf.fast5_reader.fast5_index = self.fname_mapping_filename

                self._load_index()

            elif self.prms.mode == self.WRITE_MODE:
                self.read_coords = list()

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

        self.index = load_index(self.prms.index_prefix)

        self.model = PoreModel(self.conf.pore_model)

        self.read_aln = None

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

    def init_read_aln(self, read_id):
        paf = self.mm2s[read_id]
        self.read_aln = ReadAln(self.index, paf, is_rna=self.conf.is_rna)

        if self.in_mem:
            self.read_ids = {read_id}
        else:
            self.read_ids.add(read_id)

    def get_aln_kmers(self, aln=None):
        if aln is None:
            aln = self.read_aln
        return pd.Series(
            self.index.get_kmers(aln.refmir_start, aln.refmir_end, aln.is_rna),
            pd.RangeIndex(aln.refmir_start+nt.K-1, aln.refmir_end)
        )

    def save_read(self, fast5_fname, aln=None):
        if self.prms.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")
        
        if aln is not None:
            self.read_aln = aln

        aln_i = len(self.read_coords)
        self.read_coords.append(
            (self.read_aln.ref_id, self.read_aln.ref_start, self.read_aln.ref_end, self.read_aln.read_id, aln_i)
        )

        for name in self.read_aln.dfs:
            df = getattr(self.read_aln, name)
            df = df.drop(columns=["refmir"], errors="ignore").sort_index()
            self.hdf.put("_%d/%s" % (aln_i, name), df, format="fixed")

        aln_fname = self.aln_fname(self.read_aln.read_id)

        s = "\t".join([self.read_aln.read_id, fast5_fname, "-"]) + "\n"
        self.fname_mapping_file.write(s)


    def load_read(self, read_id=None, coords=None):
        if self.read_aln is not None and read_id == self.read_aln.read_id: 
            return self.read_aln

        if read_id is None and coords is None:
            raise ValueError("read_id or coords must be specified for Track.load_read")
        elif coords is not None:
            read_id = coords.read_id
        else:
            coords = self.hdf.select("/coords", "read_id=read_id")

        group = "/_%d" % coords.i

        mm2 = self.mm2s.get(read_id, None)
        if mm2 is None: return None

        if self.prms.ref_bounds is None:
            where = None
        else:
            start = self.prms.ref_bounds.start
            end = self.prms.ref_bounds.end
            where = "index >= self.prms.ref_bounds.start & index < self.prms.ref_bounds.end"
        self.read_aln = ReadAln(self.index, mm2, is_rna=not self.conf.read_buffer.seq_fwd, ref_bounds=self.prms.ref_bounds)

        for (path, subgroups, subkeys) in self.hdf.walk(group):
            for name in subkeys:
                df = self.hdf.select(os.path.join(group, name))#, where=where)
                self.read_aln.set_df(df, name)

        for layer in self.prms.layers:
            if not layer in self.read_aln.aln.columns:
                self.read_aln.aln[layer] = self.LAYER_FNS[layer](self, self.read_aln)
                #print(self.read_aln.bcerr)
                #print(layer, self.read_aln.aln[layer].drop_duplicates())

        return self.read_aln

    def set_layers(self, layers):
        if layers is not None:
            self.prms.layers = layers
        self.layers_split = [
            tuple(l.split(".")) if "." in l else ("aln", l) 
            for l in self.prms.layers]
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

        ref_id = self.ref_id
        where = (
            "ref_id == ref_id & "
            "(ref_start < self.ref_start & ref_end > self.ref_end)" 
        )
        coords = self.hdf.select("/coords", where=where)

        #for read_id,read in self.fname_mapping.iterrows():
        for coord in coords.itertuples():
            aln = self.load_read(coords=coord)
            if aln is None: continue 

            xs = self.ref_coords.reindex(
                self.ref_coords.index.intersection(aln.aln.index)
            )

            mask_row = np.ones(self.width)#, dtype=bool)
            mask_row[xs] = False
            mask_rows.append(mask_row)

            for layer in self.prms.layers:
                row = np.zeros(self.width)
                row[xs] = aln.aln[layer]
                read_rows[layer].append(row)

            read_meta['ref_start'].append(aln.ref_start)
            read_meta['id'].append(aln.read_id)
            read_meta['fwd'].append(aln.is_fwd)

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

        self.norms = {layer: Normalize(np.min(self[layer]), np.max(self[layer])) for layer in self.prms.layers}
        for l in ["current", "dwell"]:
            layer = self[l]
            self.norms[l].vmax = min(
                layer.max(),
                np.ma.median(layer) + 2 * layer.std()
            )

        return mat

    def write_coords(self):
        df = pd.DataFrame(
                self.read_coords, 
                columns=["ref_id", "ref_start", "ref_end", "read_id", "i"]
        ).set_index(["ref_id", "ref_start", "ref_end"]).sort_index()
        self.hdf.put("coords", df, data_columns=["read_id"], format="table")

    def close(self):
        self.write_coords()
        self.hdf.close()
        self.fname_mapping_file.close()

    @property
    def hdf_filename(self):
        return os.path.join(self.prms.path, self.HDF_FNAME)

    @property
    def config_filename(self):
        return os.path.join(self.prms.path, self.CONF_FNAME)

    @property
    def fname_mapping_filename(self):
        return os.path.join(self.prms.path, self.INDEX_FNAME)
    
    @property
    def aln_dir(self):
        if self.prms.path is None: return None
        return os.path.join(self.prms.path, self.ALN_DIR)
    
    @property
    def read_count(self):
        return len(self.fname_mapping)
    
    def aln_fname(self, read_id):
        if self.prms.path is None: return None
        return os.path.join(self.aln_dir, read_id+self.ALN_SUFFIX)

    def sort(self, layer, ref):
        if isinstance(layer, str):
            layer = self.layer_idxs[layer]
        order = np.argsort(self.mat[layer,:,ref-self.ref_start])
        self.mat = self.mat[:,order,:]
        self.reads = self.reads.iloc[order]

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
        aln_a = track_a.load_read(read)
        aln_b = track_b.load_read(read)
        print (method_compare_aln(aln_a, aln_b))

def method_compare_aln(aln_a, aln_b):
    merge = aln_a.aln.join(aln_b.aln, lsuffix="_a", rsuffix="_b").dropna().set_index("refmir_a")

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
