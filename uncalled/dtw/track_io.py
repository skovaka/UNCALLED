#!/usr/bin/env python3

import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections
import time
from collections import defaultdict
import scipy

from .db import TrackSQL
from .track import AlnTrack, LAYERS
from ..index import load_index, RefCoord, str_to_coord
from ..pore_model import PoreModel
from ..fast5 import Fast5Reader, parse_read_ids
from .. import config, nt

def parse_layers(layers):
    db_layers = list() 
    fn_layers = list() 

    if isinstance(layers, str):
        layers = layers.split(",")

    ret = list()

    parsed = set()

    for layer in layers:
        spl = layer.split(".")
        if len(spl) == 2:
            group,layer = spl
        elif len(spl) == 1:
            group = "dtw"
            layer = spl[0]
        else:
            raise ValueError("Invalid layer specifier \"%s\", must contain at most one \".\"" % layer)

        if not group in LAYERS:
            raise ValueError("Invalid layer group \"%s\". Options: %s" % (group, LAYERS.keys()))

        group_layers = LAYERS[group]

        if not layer in group_layers:
            raise ValueError("Invalid layer \"%s\". Options: %s" % (group, group_layers.keys()))

        if (group, layer) in parsed:
            continue
        parsed.add((group, layer))

        yield (group, layer)
        #if group_layers[layer].fn is None:
        #    db_layers.append((group, layer))
        #else:
        #    fn_layers.append((group, layer))

    #return db_layers, fn_layers
            

class TracksParams(config.ParamGroup):
    _name = "track_io"
TracksParams._def_params(
    ("input", None, None, "Input track(s)"),
    ("output", None, None,  "Output track"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("layers", ["current", "dwell", "model_diff"], None, "Layers to load"),
    ("refstats", None, None, "Per-reference summary statistics to compute for each layer"),
    ("refstats_layers", ["current", "dwell", "model_diff"], None, "Layers to compute refstats"),
    ("read_filter", None, None, "Only load reads which overlap these coordinates"),
    ("max_reads", None, int, "Only load reads which overlap these coordinates"),
    ("index_prefix", None, str, "BWA index prefix"),
    ("load_fast5s", bool, True, "Load fast5 files"),
    ("overwrite", False, bool, "Overwrite existing databases"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("aln_chunksize", 4000, int, "Number of alignments to query for iteration"),
    ("ref_chunksize", 10000, int, "Number of reference coordinates to query for iteration"),
    ignore_toml={"input", "output", "ref_bounds", "layers"}
)

_REFSTAT_AGGS = {
    "mean" : np.mean, 
    "stdv" : np.std, 
    "var"  : np.var,
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : np.min, 
    "max"  : np.min
}

LAYER_REFSTATS = {"min", "max", "mean", "stdv", "var", "skew", "kurt"}
COMPARE_REFSTATS = {"ks"}
ALL_REFSTATS = LAYER_REFSTATS | COMPARE_REFSTATS

class RefstatSplit:
    def __init__(self, stats, track_count):
        self.layer = [s for s in stats if s in LAYER_REFSTATS]
        self.compare = [s for s in stats if s in COMPARE_REFSTATS]

        self.layer_agg = [_REFSTAT_AGGS[s] for s in self.layer]

        if len(self.layer) + len(self.compare) != len(stats):
            bad_stats = [s for s in stats if s not in ALL_REFSTATS]
            raise ValueError("Unknown stats: %s (options: %s)" % (", ".join(bad_stats), ", ".join(ALL_REFSTATS)))

        if len(self.compare) > 0 and track_count != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(self.compare))

class Tracks:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("track_io", copy_conf=True, *args, **kwargs)

        #if isinstance(self.prms.layers, str):
        #    self.prms.layers = [self.prms.layers]

        self.set_layers(self.prms.layers)

        self.dbs = dict()

        self.track_dbs = dict()

        self.aln_tracks = list()
        self.output_tracks = dict()
        self.all = self.aln_tracks

        self.prev_fast5 = dict()
        self.prev_read = dict()

        self._load_dbs(self.prms.output, True)
        self._load_dbs(self.prms.input, False)

        self.input_track_ids = [t.id for t in self.aln_tracks]
        self.output_track_ids = [t.id for _,t in self.output_tracks.items()]

        if self.prms.index_prefix is not None:
            self.index = load_index(self.prms.index_prefix)
            self._set_ref_bounds(self.prms.ref_bounds)

        if self.prms.load_fast5s:
            fast5_reads = list()
            for _,db in self.dbs.items():
                fast5_reads.append(db.get_fast5_index(self.input_track_ids))
            fast5_reads = pd.concat(fast5_reads)
            self.fast5s = self.fast5s = Fast5Reader(index=fast5_reads, conf=self.conf)

            reads = fast5_reads["read_id"]
            self.reads = pd.Index(reads[~reads.duplicated()])

            for track in self.aln_tracks:
                track.fast5s = self.fast5s

    def _set_ref_bounds(self, ref_bounds):
        if ref_bounds is not None:
            if isinstance(ref_bounds, str):
                ref_bounds = str_to_coord(ref_bounds)
            elif isinstance(ref_bounds, tuple):
                ref_bounds = RefCoord(*ref_bounds)
            self.coords = self.index.get_coord_space(ref_bounds, self.conf.is_rna)
        else:
            self.coords = None

    def set_layers(self, layers):
        self.prms.layers = layers

        self.db_layers = list()
        self.fn_layers = list()
        for group, layer in parse_layers(layers):
            if LAYERS[group][layer].fn is None:
                self.db_layers.append((group, layer))
            else:
                self.fn_layers.append((group, layer))

    def __len__(self):
        return len(self.all)

    def __getitem__(self, i):
        return self.all[i]

    def _load_dbs(self, dbs, out):
        if dbs is None:
            return

        if isinstance(dbs, str):
            dbs = [dbs]

        for db_str in dbs:
            db_file, track_names = self._db_track_split(db_str)

            if not out and not os.path.exists(db_file):
                raise OSError("Database file not found: " + db_file)

            db = self.dbs.get(db_file, None)
            if db is None:
                db = TrackSQL(db_file)
                self.dbs[db_file] = db

            if out:
                self._init_output_tracks(db, track_names)
                db.init_write()
            else:
                self._init_input_tracks(db, track_names)
    

    def _db_track_split(self, db_str):
        spl = db_str.split(":")
        if len(spl) == 1:
            filename = db_str
            track_names = None
        elif len(spl) == 2:
            filename = spl[0]
            track_names = spl[1].split(",")
        else:
            raise ValueError("Incorrect database specifier format: " + db_str)

        return os.path.abspath(filename), track_names

    def _init_track(self, db, name):
        if name in self.track_dbs:
            raise ValueError("Cannot load multiple tracks with the same name")

        self.track_dbs[name] = db

    def _init_input_tracks(self, db, track_names):
        df = db.query_track(track_names).set_index("name").reindex(track_names)

        missing = df["id"].isnull()
        if missing.any():
            bad_names = df[missing].index
            all_names = db.query_track()["name"]
            raise ValueError("alignment track not found: \"%s\" (tracks in database: \"%s\")" %
                             ("\", \"".join(bad_names), "\", \"".join(all_names)))

        for name,row in df.iterrows():
            self._init_track(db, row.name)

            conf = config.Config(toml=row.config)
            t = AlnTrack(db, row["id"], name, row["desc"], row["groups"], conf)
            self.aln_tracks.append(t)

            self.conf.load_config(conf)

    def _init_output_tracks(self, db, track_names):
        if track_names is None:
            name = os.path.splitext(os.path.basename(db.filename))[0]
            track_names = [name]

        for name in track_names:
            self._init_track(db, name)

            self.prev_fast5[name] = (None, None)
            self.prev_read[name] = None

            track = AlnTrack(db, None, name, name, "dtw", self.conf)
            self.output_tracks[name] = track

            try:
                db.init_track(track)
            except Exception as err:
                #TODO add -f option (requires delete cascade)
                if len(db.query_track(name)) > 0:
                    raise ValueError("database already contains track named \"%s\". Specify a different name, write to a different file" % name)


                raise err

    def get_fast5_reader(self):
        #TODO needs work for multiple DBs
        if self.fast5s is None:
            for db in self.dbs.values():
                fast5_index = db.get_fast5_index()
            self.fast5s = Fast5Reader(index=fast5_index, conf=self.conf)
        return self.fast5s

    def init_alignment(self, read_id, fast5, coords, group=None, layers=None, track_name=None):
        if track_name is None:
            if len(self.output_tracks) == 1:
                track_name, = self.output_tracks
            else:
                raise ValueError("Must specify track name when using multiple output tracks")

        track = self.output_tracks[track_name]
        db = track.db

        #self.prev_aln[track_name] += 1
        aln_id = db.next_aln_id()

        if fast5 == self.prev_fast5[track_name][0]:
            fast5_id = self.prev_fast5[track_name][1]
        else:
            fast5_id = db.init_fast5(fast5)
            self.prev_fast5[track_name] = (fast5, fast5_id)

        if self.prev_read[track_name] != read_id:
            db.init_read(read_id, fast5_id)
            self.prev_read[track_name] = read_id

        if layers is not None:
            if "start" in layers.columns:
                col = "start"
            elif "sample" in layers.columns:
                col = "sample"
            else:
                raise ValueError("Must initialize alignment from DataFrame with sample or start column")

            samp_start = layers[col].min()
            samp_end = layers[col].max()
        else:
            samp_start = samp_end = None

        track.alignments = pd.DataFrame({
                "id" : [aln_id],
                "track_id" : [track.id],
                "read_id" : [read_id],
                "ref_name" : [coords.ref_name],
                "ref_start" : [coords.refs.start],
                "ref_end" : [coords.refs.stop],
                "fwd" :     [coords.fwd],
                "samp_start" : [samp_start],
                "samp_end" : [samp_end],
                "tags" : [""]}).set_index("id")

        track.db.init_alignment(track.alignments)

        track.layers = None

        if layers is not None:
            track.add_layer_group(group, layers)

        track.coords = coords

        return track #aln_id

    @property
    def input_count(self):
        return len(self.aln_tracks)


    LAYER_FNS = {
        "dwell" : (lambda self,track: 
            1000 * track.layers["dtw","length"] / self.conf.read_buffer.sample_rate),
        "model_diff" : (lambda self,track: 
            track.layers["dtw","current"] - track.model[track.kmers])
    }

    def load_refs(self, ref_bounds=None, full_overlap=None, load_mat=False):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if self.coords is None:
            raise ValueError("Must set ref bounds")

        if full_overlap is None:
            full_overlap = self.prms.full_overlap

        dbfile0,db0 = list(self.dbs.items())[0]
        alignments = db0.query_alignments(self.input_track_ids, coords=self.coords, full_overlap=full_overlap)

        ids = alignments.index.to_numpy()

        layers = dict()
        for group in ["dtw"]:
            layers[group] = db0.query_layers(self.input_track_ids, self.coords, ids)
        layers = pd.concat(layers, names=["group", "layer"], axis=1)
        
        for track in self.aln_tracks:
            track_alns = alignments[alignments["track_id"] == track.id].copy()
            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i].copy()

            track.set_data(self.coords, track_alns, track_layers)
            track.calc_layers(self.fn_layers)
            
            if load_mat:
                track.load_mat()

        self.calc_refstats()

        return self.aln_tracks

    def calc_refstats(self, verbose_refs=False, cov=False):
        if self.prms.refstats is None:
            self.refstats = None
            return None

        stats = RefstatSplit(self.prms.refstats, len(self.aln_tracks))

        refstats = dict()
        grouped = [t.layers["dtw"][self.prms.refstats_layers].groupby(level="mref") for t in self.aln_tracks]

        for track,groups in zip(self.aln_tracks, grouped):
            refstats[track.name] = groups.agg(stats.layer_agg)
            if cov:
                refstats[track.name].insert(0, "cov", groups.size())

        if len(stats.compare) > 0:
            groups_a, groups_b = grouped
            mrefs_a = self.aln_tracks[0].layers.index.unique("mref")
            mrefs_b = self.aln_tracks[1].layers.index.unique("mref")
            mrefs = mrefs_a.intersection(mrefs_b)
            cmps = {l : defaultdict(list) for l in self.prms.refstats_layers}
            for mref in mrefs:
                track_a = groups_a.get_group(mref)
                track_b = groups_b.get_group(mref)
                for layer in self.prms.refstats_layers:
                    a = track_a[layer]
                    b = track_b[layer]
                    for stat in stats.compare:
                        ks = scipy.stats.stats.ks_2samp(a,b,mode="asymp")
                        cmps[layer]["stat"].append(ks.statistic)
                        cmps[layer]["pval"].append(ks.pvalue)

            refstats["ks"] = pd.concat({k : pd.DataFrame(index=mrefs, data=c) for k,c in cmps.items()}, axis=1) 

        refstats = pd.concat(refstats, axis=1, names=["track", "layer", "stat"])

        refstats.index = self.coords.mref_to_ref_index(refstats.index, multi=verbose_refs)

        self.refstats = refstats.dropna()

        return refstats

    def iter_refs(self, ref_bounds=None):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)

        dbfile0,db0 = list(self.dbs.items())[0]
        layer_iter = db0.query_layers(
            self.input_track_ids, 
            coords=self.coords, 
            order=["mref"],
            chunksize=self.prms.ref_chunksize)

        leftovers = None
        seq_coords = None
        for chunk in layer_iter:
            if leftovers is not None:
                chunk = pd.concat([leftovers, chunk])

            chunk_mrefs = chunk.index.get_level_values("mref").unique()

            if len(chunk_mrefs) == 1:
                leftovers = chunk
                continue

            if seq_coords is not None:
                coords = seq_coords.mref_intersect(chunk_mrefs[:-1])
                
            if seq_coords is None or coords is None:
                chunk_refs = self.index.mrefs_to_ref_coord(chunk_mrefs[0], chunk_mrefs[-1], not self.conf.is_rna)
                seq_refs = RefCoord(chunk_refs.name, 0, chunk_refs.ref_len, chunk_refs.fwd)
                seq_coords = self.index.get_coord_space(seq_refs, self.conf.is_rna, kmer_shift=0, load_kmers=False)
                coords = seq_coords.mref_intersect(chunk_mrefs[:-1])

            coords.kmers = self.index.mrefs_to_kmers(coords.mrefs, self.conf.is_rna)

            i = chunk_mrefs.difference(coords.mrefs)
            leftovers = chunk.loc[i]
            chunk = chunk.drop(index=i)

            layers = pd.concat({"dtw":chunk}, names=["group", "layer"], axis=1)

            aln_ids = chunk.index.unique("aln_id").to_numpy()
            alns = db0.query_alignments(self.input_track_ids, aln_id=aln_ids)

            for track in self.aln_tracks:
                track_alns = alns[alns["track_id"]==track.id].copy()
                track_layers = layers[layers.index.isin(track_alns.index, 1)].copy()
                track.set_data(coords, track_alns, track_layers)
                track.calc_layers(self.fn_layers)

            yield (coords, self.aln_tracks)

    def load_read(self, read_id, ref_bounds=None):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        else:
            self.coords = None
        
        dbfile0,db0 = list(self.dbs.items())[0]

        alns = db0.query_alignments(
            self.input_track_ids,
            read_id=read_id,
            coords=self.coords)

        self._fill_tracks(db0, alns)

        return self.aln_tracks

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if read_filter is None:
            read_filter = self.prms.read_filter
        if max_reads is None:
            max_reads = self.prms.max_reads
        
        dbfile0,db0 = list(self.dbs.items())[0]

        aln_iter = db0.query_alignments(
            self.input_track_ids,
            read_id=read_filter,
            coords=self.coords, 
            full_overlap=full_overlap, 
            order=["read_id"],
            chunksize=self.prms.aln_chunksize)

        def _iter_reads(chunk):
            for read_id, alns in chunk.groupby("read_id"):
                overlap_groups = list()
                prev = None
                for i,aln in alns.sort_values("ref_start").iterrows():
                    if prev is None or aln["ref_name"] != prev["ref_name"] or aln["ref_start"] > prev["ref_end"]:
                        overlap_groups.append([i])
                    else:
                        overlap_groups[-1].append(i)
                    prev = aln
                    
                for group in overlap_groups:
                    self._fill_tracks(db0, alns.loc[group])
                    yield (read_id, self.aln_tracks)

        n = 0
        leftovers = pd.DataFrame()
        for chunk in aln_iter:
            if leftovers is not None:
                chunk = pd.concat([leftovers, chunk])

            end = chunk["read_id"] == chunk["read_id"].iloc[-1]
            leftovers = chunk[end]
            chunk = chunk[~end]

            for read in _iter_reads(chunk):
                n += 1
                if max_reads is not None and n >= max_reads:
                    return read
                yield read

        for read in _iter_reads(leftovers):
            n += 1
            if max_reads is not None and n >= max_reads:
                return read
            yield read

    def _fill_tracks(self, db, alns):
        ids = list(alns.index)

        layers = dict()
        #TODO parse which track layers to import
        for group in ["dtw"]: #self.groups:
            layers[group] = db.query_layers(self.input_track_ids, self.coords, ids) #, index=["aln_id","mref"])

        layers = pd.concat(layers, names=["group", "layer"], axis=1)

        for track in self.aln_tracks:
            track_alns = alns[alns["track_id"] == track.id]

            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i].copy()

            name = track_alns["ref_name"].iloc[0]
            fwd = track_alns["fwd"].iloc[0]
            start = track_alns["ref_start"].min()
            end = track_alns["ref_end"].max()
            ref_coord = RefCoord(name, start, end, fwd)
            track_coords = self.index.get_coord_space(
                ref_coord, self.conf.is_rna, kmer_shift=0, load_kmers=True)

            track.set_data(track_coords, track_alns, track_layers)
            track.calc_layers(self.fn_layers)
    
    def close(self):
        for filename, db in self.dbs.items():
            db.close()
