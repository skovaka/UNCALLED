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
    ignore_toml={"mode", "overwrite"}
)

class AlnTrack:
    DB_FNAME = "alns.db"
    CONF_FNAME = "conf.toml"
    ALN_DIR = "alns"

    INDEX_FNAME = "filename_mapping.txt"
    INDEX_HEADER = "read_id\tfilename\taln_file"

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
        self.aln_reads = dict()

        if self.prms.name is None:
            self.prms.name = os.path.basename(self.prms.path)

        if not self.in_mem:
            if self.prms.mode == self.WRITE_MODE:
                os.makedirs(self.prms.path, exist_ok=self.prms.overwrite)

            self.fname_mapping_file = open(self.fname_mapping_filename, self.prms.mode)

            self.con = sqlite3.connect(self.db_filename)

            if self.prms.mode == self.READ_MODE:
                self.conf.load_toml(self.config_filename)

                if len(self.conf.fast5_reader.fast5_index) == 0:
                    self.conf.fast5_reader.fast5_index = self.fname_mapping_filename

                self._load_fast5_index()

        self.index = load_index(self.prms.index_prefix)
        self.model = PoreModel(self.conf.pore_model)
        self.read_aln = None
        self.mat = None

        self.set_ref_bounds(self.prms.ref_bounds)

        if self.conf.dtw.mm2_paf is not None:
            read_filter = set(self.conf.fast5_reader.read_filter)

        if not self.in_mem:
            if self.prms.load_mat and self.prms.ref_bounds is not None:
                self.load_region(self.prms.ref_bounds)

            elif self.prms.mode == self.WRITE_MODE:
                self.read_coords = list()

                self.conf.to_toml(self.config_filename)
                self.fname_mapping_file.write(self.INDEX_HEADER + "\n")
                self._init_tables()

    def _init_tables(self):
        cur = self.con.cursor()
        cur.execute("""
            CREATE TABLE IF NOT EXISTS "%s.alns" (
                aln_id INTEGER PRIMARY KEY,
                read_id TEXT,
                ref_name TEXT,
                ref_start INTEGER,
                ref_end INTEGER,
                fwd INTEGER
            );""" % self.prms.name)
        cur.execute("""
            CREATE TABLE IF NOT EXISTS "%s.dtw" (
                mref INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                current REAL,
                PRIMARY KEY (mref, aln_id),
                FOREIGN KEY (aln_id) REFERENCES "%s.dtw" (aln_id)
            );""" % (self.prms.name, self.prms.name))
        self.con.commit()

    def _init_aln(self, aln):
        cur = self.con.cursor()
        cur.execute("INSERT INTO \"%s.alns\" VALUES (?,?,?,?,?,?)" % self.prms.name,
                     (aln.aln_id, aln.read_id, aln.ref_name, aln.ref_start, aln.ref_end, aln.is_fwd))
        print(cur.lastrowid)
                     #("\"%s.alns\"" % self.prms.name, aln.id, aln.read_id, aln.ref_name, aln.ref_start, aln.ref_end, aln.is_fwd)))

    def __contains__(self, read_id):
        return read_id in self.read_ids

    def _load_fast5_index(self):
        self.fname_mapping = pd.read_csv(self.fname_mapping_file, sep="\t", index_col="read_id")
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
            self.ref_coords = None
            self.mref_coords = None
            return

        ref_index = pd.RangeIndex(self.ref_start, self.ref_end-nt.K+1)

        self.width = len(ref_index)
        self.height = None

        ref_coords = list()
        kmers = list()

        for fwd in [False, True]:
            mref = self._ref_coords_to_mrefs(ref_bounds, fwd)
            kmer = self.index.get_kmers(mref.start, mref.stop+nt.K-1, fwd)

            if fwd == self.conf.is_rna:
                mref = mref[::-1]
                kmer = kmer[::-1]

            ref_coords.append(
                pd.DataFrame(index=mref, data={"ref" : ref_index, "fwd" : fwd})
            )
            kmers.append(pd.Series(index=mref, data=kmer))

        #TODO rename coords to something better
        self.ref_coords = pd.concat(ref_coords).sort_index()
        self.mref_coords = self.ref_coords.reset_index().set_index(["fwd", "ref"])["mref"].unstack(level=0)
        self.kmers = pd.concat(kmers).sort_index()

        print(self.mref_coords)

        #print("MREF")
        #print(self.mrefs)
        #print("COORD")
        #print(self.coords)

            

    def init_read_aln(self, read_id, bounds):
        aln_id = len(self.aln_reads)

        if isinstance(bounds, RefCoord):
            bounds = self._ref_coords_to_mrefs(bounds)
        elif not isinstance(bounds, pd.RangeIndex):
            raise ValueError("ReadAlns can only be initialized with RangeIndex or RefCoord bounds")

        if self.ref_coords is not None:
            fwd = self.index.is_mref_fwd(bounds.min(), self.conf.is_rna)
            bounds = bounds.intersection(self._bounds[fwd])
            if len(bounds) == 0: return False

        self.read_aln = ReadAln(aln_id, read_id, bounds, index=self.index, is_rna=self.conf.is_rna)
        self._init_aln(self.read_aln)

        self.aln_reads[aln_id] = read_id

        if self.in_mem:
            self.read_ids = {read_id}
        else:
            self.read_ids.add(read_id)

        return True

    def load_aln_kmers(self, aln=None, store=True):
        if aln is None:
            aln = self.read_aln

        kmers = pd.Series(
            self.index.get_kmers(aln.mref_start-nt.K+1, aln.mref_end, aln.is_rna),
            #aln.mrefs
            pd.RangeIndex(aln.mref_start, aln.mref_end)
        )

        if store:
            self.read_aln.aln["kmer"] = kmers
        return kmers

    def save_read(self, fast5_fname, aln=None):
        if self.prms.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")
        
        if aln is not None:
            self.read_aln = aln

        #columns=["aln_id", "read_id", "ref_id", "ref_start", "ref_end", "fwd", "primary"]
        aln_i = len(self.read_coords)
        self.read_coords.append(
            (aln_i, self.read_aln.read_id, self.read_aln.ref_id, self.read_aln.ref_start, self.read_aln.ref_end, self.read_aln.is_fwd, True)
        )

        #for name in self.read_aln.dfs:
        #    df = getattr(self.read_aln, name)
        #    df = df.drop(columns=["mref"], errors="ignore").sort_index()
        #    self.hdf.put("_%d/%s" % (aln_i, name), df, format="fixed")
        df = self.read_aln.aln
        #print(df)
        df.to_sql("%s.dtw" % self.prms.name, self.con, if_exists="append", index=True, index_label="mref")

        s = "\t".join([self.read_aln.read_id, fast5_fname, "-"]) + "\n"
        self.fname_mapping_file.write(s)


    def load_read(self, read_id=None, coords=None, load_kmers=True):
        if self.read_aln is not None and read_id == self.read_aln.read_id: 
            return self.read_aln

        if read_id is None and coords is None:
            raise ValueError("read_id or coords must be specified for AlnTrack.load_read")
        elif coords is not None:
            read_id = coords.read_id
        else:
            coords = self.hdf.select("/coords", "read_id=read_id").iloc[0]

        aln_id = coords.aln_id
        group = "/_%d" % aln_id

        #if self.prms.ref_bounds is None:
        #    where = None
        #else:
        #    start = self.prms.ref_bounds.start
        #    end = self.prms.ref_bounds.end
        #    where = "index >= self.prms.ref_bounds.start & index < self.prms.ref_bounds.end"

        if self.mref_coords != None:
            mrefs = self.mref_coords[coords.fwd]
            if len(mrefs) == 0: return False
        else:
            mrefs = None

        self.read_aln = ReadAln(aln_id, read_id, mrefs, index=self.index, is_rna=self.conf.is_rna)

        for (path, subgroups, subkeys) in self.hdf.walk(group):
            for name in subkeys:
                df = self.hdf.select(os.path.join(group, name))#, where=where)
                self.read_aln.set_df(df, name)

        if not self.read_aln.empty:
            if load_kmers:
                self.load_aln_kmers()

            for layer in self.prms.layers:
                if not layer in self.read_aln.aln.columns:
                    self.read_aln.aln[layer] = self.LAYER_FNS[layer](self, self.read_aln)

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


        t = time.time()

        select = "SELECT mref, \"%s.dtw\".aln_id, start, length, current FROM \"%s.dtw\"" 
        where = " WHERE (mref >= ? AND mref <= ?)"
        params = [int(self.mref_coords[True].min()), str(self.mref_coords[True].max())]
        name_count = 2

        if self.prms.full_overlap:
            select = select + " JOIN \"%s.alns\" ON \"%s.alns\".aln_id = \"%s.dtw\".aln_id"
            where = where + " AND (ref_start < ? AND ref_end > ?)"
            params += [int(self.mref_coords.index[0]), int(self.mref_coords.index[-1])]
            name_count += 3


        query = select + where
        fmt = ((self.prms.name,)*name_count)

        self.df = pd.read_sql_query(
            query % fmt, self.con, params=params,
        )#.rename(index=self.mref_to_ref) \
         #.rename_axis("ref", axis=0).sort_index()

        print(query)
        print(params)
        print(self.df)

        ids = self.df["aln_id"][~self.df["aln_id"].duplicated()].to_numpy(dtype=str)

        #TODO deal with multiple strands
        #add OR clause to WHERE, compute ref, pivot, return one df/mat?
        #or two queries, pivot seperately, return two dfs/mats?
        self.reads = pd.read_sql_query(
            "SELECT * FROM \"%s.alns\" WHERE aln_id IN (%s)" 
            % (self.prms.name, ",".join(["?"]*len(ids))),
            self.con, params=ids#list(self.mat.index)
        ).rename(columns={"read_id" : "id"}).sort_values("ref_start")

        self.df.set_index(["mref", "aln_id"], inplace=True)

        self.mat = self.df.reset_index().pivot(index="aln_id", columns="mref") \
                   .rename(columns=self.ref_coords["ref"]) \
                   .rename_axis(("layer","ref"), axis=1) \
                   .sort_index(axis=1,level=1) 


        self.mat = self.mat.reindex(self.reads["aln_id"], copy=False)

        self.has_fwd = np.any(self.reads['fwd'])
        self.has_rev = not np.all(self.reads['fwd'])

        self.height = len(self.reads)


        return self.mat

    def add_layer(self, name, mat):
        self.mat = np.ma.concatenate([self.mat, [mat]])
        self.set_layers(self.prms.layers + [name])

    def get_pileup(self, layer):
        return np.flip(np.sort(self[layer], axis=0), axis=0)

    def write_coords(self):
        df = pd.DataFrame(
                self.read_coords, 
                columns=["aln_id", "read_id", "ref_id", "ref_start", "ref_end", "fwd", "primary"]
        ).set_index(["ref_id", "ref_start", "ref_end"]).sort_index()
        self.hdf.put("coords", df, data_columns=["read_id"], format="table")

    def close(self):
        if self.in_mem: return
        #if self.prms.mode == "w":
        #    self.write_coords()
        #self.hdf.close()
        self.con.close()
        self.fname_mapping_file.close()

    #@property
    #def name(self):
    #    if self.prms.path is None:
    #        return self.fast5s
    #    return self.prms.path.split("/")[-1]

    @property
    def db_filename(self):
        return os.path.join(self.prms.path, self.DB_FNAME)

    @property
    def config_filename(self):
        return os.path.join(self.prms.path, self.CONF_FNAME)

    @property
    def fname_mapping_filename(self):
        return os.path.join(self.prms.path, self.INDEX_FNAME)
    
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

        for i,l in enumerate(AlnTrack.CMP_LAYERS):
            for j,rf in enumerate(self.mref_coords[True]):
                a = self[l,:,rf]
                b = track_b[l,:,rf]
                ks = scipy.stats.mstats.ks_2samp(a,b,mode="asymp")
                ks_stats[i][j] = ks[0]

        return ks_stats

    def calc_pca(self, layer, ref_start, ref_end, n_components=2):
        x = self[layer,:,ref_start:ref_end].T
        pc = PCA(n_components=n_components).fit_transform(x)
        data = {"read_id" : self.reads["id"]}
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
        print (method_compare_aln(aln_a, aln_b))

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
