#!/usr/bin/env python3

import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections
import time
from collections import defaultdict
import scipy

from .db import TrackSQL, delete, Eventalign
from .aln_track import AlnTrack, LAYERS, parse_layers
from ..index import load_index, RefCoord, str_to_coord
from ..pore_model import PoreModel
from ..fast5 import Fast5Reader, parse_read_ids
from .. import config, nt

            
class TracksParams(config.ParamGroup):
    _name = "tracks"
TracksParams._def_params(
    ("input", None, None, "Input tracks specifier. Should be in the format <file.db>[:<track1>[,<track2>...]]. If no track names are specified, all tracks will be loaded from the database."),
    ("output", None, None,  "Output track specifier. Should be in the format <file.db>[:<track_name>], where file.db is the output sqlite database. If <track_name> is not specified, will be same as filename (without extension)"),
    ("output_format", "db", str,  "Output format (db, nanopolish)"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("layers", ["dtw","bcaln","cmp","bc_cmp"], None, "Layers to load (e.g. current, dwell, model_diff)"),
    ("refstats", None, None, "Per-reference summary statistics to compute for each layer"),
    ("refstats_layers", None, None, "Layers to compute refstats"),
    ("read_filter", None, None, "Only load reads which overlap these coordinates"),
    ("max_reads", None, int, "Only load reads which overlap these coordinates"),
    ("index_prefix", None, str, "BWA index prefix"),
    ("load_fast5s", bool, True, "Load fast5 files"),
    ("overwrite", False, bool, "Overwrite existing tracks"),
    ("append", False, bool, "Append reads to existing tracks"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("min_coverage", 1, int, "Reference positions with less than this coverage will be excluded from each track (or all tracks if shared_refs_only is true)"),
    ("shared_reads_only", False, bool, "If true will only contain reads shared between all tracks"),
    ("shared_refs_only", False, bool, "If true will only contain reference positions where all tracks have sufficient coverage (see min_coverage)"),
    ("load_mat", False, bool, "If true will pivot layers into a matrix"), #TODO change to mat_layers, only do it for them
    ("aln_chunksize", 4000, int, "Number of alignments to query for iteration"),
    ("ref_chunksize", 10000, int, "Number of reference coordinates to query for iteration"),
    ignore_toml={"input", "output", "output_format", "ref_bounds", "layers", "full_overlap", "overwrite", "append","refstats", "refstats_layers", "read_filter"}
)

_REFSTAT_AGGS = {
    "mean" : np.mean, 
    "median" : np.median, 
    "q5" : (lambda x: np.quantile(x, 0.05)),
    "q95" : (lambda x: np.quantile(x, 0.95)),
    "q25" : (lambda x: np.quantile(x, 0.25)),
    "q75" : (lambda x: np.quantile(x, 0.75)),
    "stdv" : np.std, 
    "var"  : np.var,
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : np.min, 
    "max"  : np.min,
}

LAYER_REFSTATS = {"min", "max", "mean", "median", "stdv", "var", "skew", "kurt", "q25", "q75", "q5", "q95"}
COMPARE_REFSTATS = {"ks"}
ALL_REFSTATS = LAYER_REFSTATS | COMPARE_REFSTATS

REFSTAT_LABELS = {
    "min" : "Minimum", 
    "max" : "Maximum", 
    "mean" : "Mean", 
    "median" : "Median", 
    "stdv" : "Std. Dev.", 
    "var" : "Variance", 
    "skew" : "Skewness", 
    "kurt" : "Kurtosis",
    "ks" : "KS",
    "q5" : "5% Quantile",
    "q25" : "25% Quantile",
    "q75" : "75% Quantile",
    "q95" : "95% Quantile",
}

BUILTIN_TRACKS = {"_refstats", "_dtwstats", "_readstats"}

CMP_GROUPS = {"cmp", "bc_cmp"}

class RefstatsSplit:
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
        if len(args) > 0 and isinstance(args[0], Tracks):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

    def _init_new(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("tracks", copy_conf=True, *args, **kwargs)

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers, add_deps=False))

        #TODO refactor into TrackSQL, abstract into self.io
        self.dbs = dict()
        self.track_dbs = dict() #TODO do we need this?

        self.alns = list()
        self.output_tracks = dict()
        self._tracks = dict()

        self.prev_fast5 = dict()
        self.prev_read = dict()

        if self.prms.output_format == "db":
            self._load_dbs(self.prms.output, True)
            self._load_dbs(self.prms.input, False)
        else:
            self.io = Eventalign(self.prms.output, "wb")

            name = ""
            track = AlnTrack(None, None, name, name, self.conf)
            self.output_tracks[name] = track

        if self.prms.index_prefix is None:
            raise RuntimeError("No reference index")

        self._aln_track_ids = [t.id for t in self.alns]

        self.index = load_index(self.prms.index_prefix)

        self.coords = self._ref_bounds_to_coords(self.prms.ref_bounds)

        if self.coords is not None:
            self.load()

        if self.prms.load_fast5s and len(self.dbs) > 0:
            fast5_reads = list()
            for _,db in self.dbs.items():
                fast5_reads.append(db.get_fast5_index(self._aln_track_ids))
            fast5_reads = pd.concat(fast5_reads)
            files = fast5_reads["filename"].unique()
            self.fast5s = Fast5Reader(
                index=fast5_reads, 
                conf=self.conf)

            reads = fast5_reads["read_id"]
            self.reads = pd.Index(reads[~reads.duplicated()])

            for track in self.alns:
                track.fast5s = self.fast5s

    def _init_slice(self, parent, coords, tracks):
        self.conf = parent.conf 
        self.prms = parent.prms

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers, add_deps=False))

        self.dbs = parent.dbs
        self.track_dbs = parent.track_dbs
        self.output_tracks = parent.output_tracks
        self.prev_fast5 = parent.prev_fast5
        self.prev_read = parent.prev_read
        self.index = parent.index
        self.fast5s = parent.fast5s
        self._aln_track_ids = parent._aln_track_ids
        self.refstats = None

        self.coords = coords
        self._tracks = tracks

        self.alns = list()
        if self._tracks is not None:
            for name,track in tracks.items():
                if name in BUILTIN_TRACKS:
                    setattr(self, name[1:], track)
                elif isinstance(track, AlnTrack):
                    self.alns.append(track)

    @property
    def all_empty(self):
        for t in self.alns:
            if not t.empty:
                return False
        return True

    @property
    def any_empty(self):
        for t in self.alns:
            if t.empty:
                return True
        return False

    def _ref_bounds_to_coords(self, ref_bounds):
        if ref_bounds is not None:
            if isinstance(ref_bounds, str):
                ref_bounds = str_to_coord(ref_bounds)
            elif isinstance(ref_bounds, tuple):
                ref_bounds = RefCoord(*ref_bounds)
            return self.index.get_coord_space(ref_bounds, self.conf.is_rna)
        return None

    def set_layers(self, layers):
        self.prms.layers = layers

        self.db_layers = list()
        self.fn_layers = list()
        self.cmp_layers = list()
        for group, layer in parse_layers(layers):
            if LAYERS[group][layer].fn is None:
                if group in CMP_GROUPS:
                    self.cmp_layers.append((group, layer))
                else:
                    self.db_layers.append((group, layer))
            else:
                self.fn_layers.append((group, layer))

    def __len__(self):
        return len(self.alns)

    def __getitem__(self, i):
        return self.alns[i]

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
            raise ValueError("Invalid database specifier format: " + db_str)

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
            t = AlnTrack(db, row["id"], name, row["desc"], conf)
            self.alns.append(t)
            self._tracks[t.name] = t

            self.conf.load_config(conf)

    def _init_output_tracks(self, db, track_names):
        if track_names is None:
            name = os.path.splitext(os.path.basename(db.filename))[0]
            track_names = [name]

        for name in track_names:
            self._init_track(db, name)

            self.prev_fast5[name] = (None, None)
            self.prev_read[name] = None

            track = AlnTrack(db, None, name, name, self.conf)
            self.output_tracks[name] = track

            try:
                db.init_track(track)
            except Exception as err:
                if len(db.query_track(name)) > 0:
                    if self.prms.append:
                        continue
                    elif self.prms.overwrite:
                        sys.stderr.write("Deleting existing track...\n")
                        delete(name, db)
                        db.init_track(track)
                        continue
                    else:
                        raise ValueError("database already contains track named \"%s\". Specify a different name, write to a different file" % name)


                raise err

    def aln_layers(self, layer_filter=None):
        ret = pd.Index([])
        for track in self.alns:
            layers = track.layers.columns
            if layer_filter is not None:
                layers = layers.intersection(layer_filter)
            ret = ret.union(layers)
        return ret
            

    
    def get_fast5_reader(self):
        if self.fast5s is None:
            for db in self.dbs.values():
                fast5_index = db.get_fast5_index()
            self.fast5s = Fast5Reader(index=fast5_index, conf=self.conf)
        return self.fast5s

    def _track_or_default(self, track_name):
        if track_name is None:
            if len(self.output_tracks) == 1:
                track_name, = self.output_tracks
            else:
                raise ValueError("Must specify track name when using multiple output tracks")
        return self.output_tracks[track_name]

    def collapse_events(self, dtw):
        dtw["cuml_mean"] = dtw["length"] * dtw["current"]

        grp = dtw.groupby("mref")

        mrefs = grp["mref"].first()

        lengths = grp["length"].sum()

        dtw = pd.DataFrame({
            "mref"    : mrefs.astype("int64"),
            "start"  : grp["start"].min().astype("uint32"),
            "length" : lengths.astype("uint32"),
            "current"   : grp["cuml_mean"].sum() / lengths
        })

        return dtw.set_index("mref").sort_index()

    def write_events(self, events, track_name=None, aln_id=None):
        if self.prms.output_format == "db":
            dtw = self.collapse_events(events)
            self.write_layers({"dtw" : dtw}, track_name, aln_id)
        elif self.prms.output_format == "eventalign":
            track = self._track_or_default(track_name)
            self.io.write_events(track, events)
        
    def write_layers(self, layers, track_name=None, aln_id=None):
        track = self._track_or_default(track_name)
        df = track.add_layer_group(layers, aln_id)

        if self.prms.output_format == "db":
            db = self.track_dbs[track.name]
            db.write_layers(df)

    def write_alignment(self, read_id, fast5, coords, layers={}, track_name=None):
        track = self._track_or_default(track_name)

        if len(self.track_dbs) > 0:
            db = self.track_dbs[track.name]
            if fast5 == self.prev_fast5[track.name][0]:
                fast5_id = self.prev_fast5[track.name][1]
            else:
                fast5_id = db.init_fast5(fast5)
                self.prev_fast5[track.name] = (fast5, fast5_id)

            if self.prev_read[track.name] != read_id:
                db.init_read(read_id, fast5_id)
                self.prev_read[track.name] = read_id
        else:
            db = self.io

        aln_id = db.next_aln_id()


        if len(layers) > 0:
            starts = None
            for group,vals in layers.items():
                if "start" in vals.columns:
                    starts = vals["start"]
                    break

            if starts is None:
                raise ValueError("Must initialize alignment from DataFrame with start column")

            samp_start = starts.min()
            samp_end = starts.max()
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
        
        db.init_alignment(track.alignments)
        #aln_id = db.next_aln_id()
        #fast5_id = db.init_fast5(fast5)
        #db.init_read(read_id, fast5_id)

        track.layers = None

        track.coords = coords

        for group,vals in layers.items():
            if group == "dtw":
                self.write_events(vals, track.name, aln_id)
            else:
                self.write_layers({group : vals}, track.name, aln_id)

        return aln_id, coords

    @property
    def input_count(self):
        return len(self.alns)

    def _verify_read(self):
        if len(self._aln_track_ids) == 0:
            raise RuntimeError("No input tracks have been loaded")

    #TODO read_ids, track_alns, max_cov(?)
    def slice(self, ref_start=None, ref_end=None, reads=None, order=["fwd","ref_start"]):
        if self.coords is None:
            raise IndexError("Cannot slice empty Tracks")

        hasbounds = (ref_start is not None, ref_end is not None)
        if np.all(hasbounds):
            coords = self.coords.ref_slice(ref_start, ref_end)
        elif np.any(hasbounds):
            raise IndexError(f"Invalid bounds {ref_start}-{ref_end}")
        else:
            coords = self.coords

        tracks = dict()
        for name,track in self._tracks.items():
            if isinstance(track, pd.DataFrame):
                tracks[name] = track.loc[coords.refs]

            elif isinstance(track, AlnTrack):
                tracks[name] = track.slice(coords, reads=reads, order=order)

        return Tracks(self, coords, tracks)

    def get_shared_reads(self):
        read_ids = pd.Index(self.alns[0].read_ids)
        for track in self.alns[1:]:
            read_ids = read_ids.intersection(track.read_ids)
        return read_ids

    def get_all_reads(self):
        read_ids = pd.Index(self.alns[0].read_ids)
        for track in self.alns[1:]:
            if not track.empty:
                read_ids = read_ids.union(track.read_ids)
        return read_ids

    def slice_shared_reads(self):
        return self.slice(reads=self.get_shared_reads(), order="read_id")
            
    def load(self, ref_bounds=None, full_overlap=None, read_filter=None, load_mat=False):
        self._verify_read()

        if ref_bounds is not None:
            self.coords = self._ref_bounds_to_coords(ref_bounds)

        if read_filter is None:
            read_filter = self.prms.read_filter

        if self.coords is None:
            raise ValueError("Must set ref bounds")

        if full_overlap is None:
            full_overlap = self.prms.full_overlap

        if load_mat is None:
            load_mat = self.prms.load_mat

        dbfile0,db0 = list(self.dbs.items())[0]

        layers = db0.query_layers(self.db_layers, self._aln_track_ids, self.coords, full_overlap=full_overlap, read_id=read_filter)

        ids = layers.index.get_level_values("aln_id").unique().to_numpy()

        alignments = db0.query_alignments(aln_id=ids)#self._aln_track_ids, coords=self.coords, full_overlap=full_overlap)


        for track in self.alns:
            track_alns = alignments[alignments["track_id"] == track.id]
            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i]

            track.set_data(self.coords, track_alns, track_layers)
            track.calc_layers(self.fn_layers)
        
        self.load_compare(ids)
        self.calc_refstats()

        if load_mat:
            for track in self.alns:
                track.load_mat()

        return self.alns

    def load_compare(self, aln_ids=None):
        if len(self.cmp_layers) == 0:
            return

        dbfile0,db = list(self.dbs.items())[0]

        self.cmp = db.query_compare(self.cmp_layers, self._aln_track_ids, self.coords, aln_ids)
        #TODO add aln_ids

        groups = self.cmp.index.get_level_values("group_b").unique()
        if "bcaln" in groups:
            bcalns = self.cmp.loc[(slice(None), slice(None), slice(None), "bcaln"),:]
        else:
            bcalns = None
        if "dtw" in groups:
            dtws = self.cmp.loc[(slice(None), slice(None), slice(None), "dtw"),:]
        else:
            dtws = None

        for track in self.alns:
            def _add_group(group, df):
                df = df.reset_index(["aln_b", "group_b"])
                df = df[df.index.get_level_values("mref").isin(track.layer_mrefs)]
                df.rename(index=track.coords.mref_to_ref, level=0, inplace=True)
                df.index.names = ["ref", "aln_id"]
                #print("A", list(track.layers.columns))
                df = pd.concat({group : df.reindex(track.layers.index)}, axis=1)
                track.layers = pd.concat([track.layers, df], axis=1).dropna(axis=1,how="all")

            try:
                if bcalns is not None:
                    _add_group("bc_cmp", bcalns)
                if dtws is not None:
                    _add_group("cmp", dtws)
            except:
                sys.stderr.write("Failed to write compare group\n")
                sys.stderr.write(str(track.alignments))

    def calc_compare(self, group_b, calc_jaccard, calc_mean_ref_dist, save):
        if len(self.alns) == 0:
            raise ValueError("Must input at least one track")

        cols = self.alns[0].layers.columns.get_level_values("group").unique()
        if (group_b == "dtw" and "cmp" in cols) or (group_b == "bcaln" and "bc_cmp" in cols):
            sys.stderr.write(f"Read already has compare group. Skipping\n")
            return None

        if group_b == "dtw":
            if len(self.alns) != 2:
                raise ValueError("Must input exactly two tracks to compare dtw alignments")

            df = self.alns[0].cmp(self.alns[1], calc_jaccard, calc_mean_ref_dist)

        elif group_b == "bcaln":
            if len(self.alns) > 2:
                raise ValueError("Must input one or two tracks to compare dtw to bcaln")
            
            if len(self.alns) == 2:
                df = self.alns[0].bc_cmp(self.alns[1])
            else:
                df = self.alns[0].bc_cmp()

        df = df.dropna(how="all")

        if save:
            df.rename(index=lambda r: self.alns[0].coords.ref_to_mref(r, True), level=0, inplace=True)
            self.alns[0].db.write_layers(
                pd.concat({"cmp" : df}, names=["group", "layer"], axis=1), 
                index=["mref", "aln_a", "aln_b", "group_b"])
        else:
            print(df.to_csv(sep="\t"))

    def calc_refstats(self, verbose_refs=False, cov=False):
        if self.prms.refstats is None or len(self.prms.refstats) == 0 or len(self.prms.refstats_layers) == 0 or self.all_empty:
            self.refstats = None
            return None

        stats = RefstatsSplit(self.prms.refstats, len(self.alns))

        refstats = dict()
        grouped = [
            t.layers[self.prms.refstats_layers].groupby(level=0)
            for t in self.alns]

        for track,groups in zip(self.alns, grouped):
            if track.empty:
                refstats[track.name] = None
                continue


            refstats[track.name] = groups.agg(stats.layer_agg)
            rename = ({
                old[-1] : new
                for old,new in zip(refstats[track.name].columns, stats.layer)
            })
            refstats[track.name].rename(columns=rename, inplace=True)

            if cov:
                refstats[track.name].insert(0, "cov", groups.size())

        if len(stats.compare) > 0:
            if self.any_empty:
                refstats["ks"] = None
            
            else:
                groups_a, groups_b = grouped
                refs_a = self.alns[0].layers.index.unique("ref")
                refs_b = self.alns[1].layers.index.unique("ref")
                refs = refs_a.intersection(refs_b)
                cmps = {l : defaultdict(list) for l in self.prms.refstats_layers}
                for ref in refs:
                    track_a = groups_a.get_group(ref)
                    track_b = groups_b.get_group(ref)
                    for layer in self.prms.refstats_layers:
                        a = track_a[layer]
                        b = track_b[layer]
                        for stat in stats.compare:
                            ks = scipy.stats.stats.ks_2samp(a,b,mode="asymp")
                            cmps[layer]["stat"].append(ks.statistic)
                            cmps[layer]["pval"].append(ks.pvalue)

                refstats["ks"] = pd.concat({k : pd.DataFrame(index=refs, data=c) for k,c in cmps.items()}, axis=1) 
        
        if np.any([df is None for df in refstats.values()]):
            columns = None
            for df in refstats.values():
                if df is not None:
                    columns = df.columns
                    break
            for name,df in refstats.items():
                if df is None:
                    refstats[name] = pd.DataFrame(columns=columns)
        refstats = pd.concat(refstats, axis=1, names=["track", "group", "layer", "stat"])

        #TODO make verbose ref indexing
        #refstats.index = self.alns[0].coords.mref_to_ref_index(refstats.index, multi=verbose_refs)

        self.refstats = refstats.dropna()

        self._tracks["_refstats"] = self.refstats

        return refstats

    def iter_refs(self, ref_bounds=None):
        if ref_bounds is not None:
            coords = self._ref_bounds_to_coords(ref_bounds)
        else:
            coords = self.coords

        t0 = time.time()
        dbfile0,db0 = list(self.dbs.items())[0]
        layer_iter = db0.query_layers(
            self.db_layers, 
            self._aln_track_ids, 
            coords=coords, 
            order=["mref"],
            chunksize=self.prms.ref_chunksize)

        t0 = time.time()

        def get_full_coords(mref):
            chunk_refs = self.index.mrefs_to_ref_coord(mref, mref, not self.conf.is_rna)
            seq_refs = RefCoord(chunk_refs.name, 0, chunk_refs.ref_len, chunk_refs.fwd)
            return self.index.get_coord_space(seq_refs, self.conf.is_rna, load_kmers=False)

        def next_coords(seq_coords, mrefs):
            if len(mrefs) > 1:
                mrefs = mrefs[:-1]

            coords = seq_coords.mref_intersect(mrefs)
            if coords is None:
                seq_coords = get_full_coords(mrefs[0])
                coords = seq_coords.mref_intersect(mrefs)

            return seq_coords, coords

        chunk = next(layer_iter)
        seq_coords = get_full_coords(chunk.index.get_level_values("mref")[0])

        t0 = time.time()
                
        while len(chunk) > 0:

            if len(chunk.index.get_level_values("mref").unique()) <= 1:
                chunk = pd.concat([chunk, next(layer_iter, pd.DataFrame())])

            chunk_mrefs = chunk.index.get_level_values("mref").unique()

            seq_coords, coords = next_coords(seq_coords, chunk_mrefs)

            coords.set_kmers(self.index.mrefs_to_kmers(coords.mrefs, self.conf.is_rna, False))

            i = chunk_mrefs.difference(coords.mrefs)
            leftovers = chunk.loc[i]
            layers = chunk.drop(index=i)

            aln_ids = layers.index.unique("aln_id").to_numpy()
            alns = db0.query_alignments(self._aln_track_ids, aln_id=aln_ids)

            ret = self._tables_to_tracks(coords, alns, layers)

            chunk = leftovers

            if not ret.all_empty:
                yield ret

    def _mref_to_coords(self, mref):
        pass

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None):
        
        if ref_bounds is not None and not isinstance(ref_bounds, RefCoord):
            ref_bounds = RefCoord(ref_bounds)

        if (self.coords is None or
            (read_filter is not None and 
             len(self.get_all_reads().intersection(read_filter)) < len(read_filter)) or
            (ref_bounds is not None and not self.coords.contains(ref_bounds))):
            gen = self.iter_reads_db(read_filter, ref_bounds, full_overlap, max_reads)
        else:
            gen = self.iter_reads_slice(read_filter, ref_bounds)

        for read_id,chunk in gen:
            yield read_id,chunk
            
    def iter_reads_slice(self, reads=None, ref_bounds=None):
        all_reads = self.get_all_reads()
        if reads is not None:
            all_reads = all_reads.intersection(reads)

        if ref_bounds is None:
            ref_start = ref_end = None
        else:
            ref_start = ref_bounds.start
            ref_end = ref_bounds.end
        
        for read_id in all_reads:
            yield read_id, self.slice(ref_start, ref_end, [read_id])

    def iter_reads_db(self, reads, ref_bounds, full_overlap, max_reads):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if reads is None:
            reads = self.prms.read_filter
        if max_reads is None:
            max_reads = self.prms.max_reads
        
        dbfile0,db0 = list(self.dbs.items())[0]

        layer_iter = db0.query_layers(
            self.db_layers, 
            self._aln_track_ids,
            read_id=reads,
            coords=self.coords, 
            full_overlap=full_overlap, 
            order=["read_id"],
            chunksize=self.prms.ref_chunksize)

        aln_leftovers = pd.DataFrame()
        layer_leftovers = pd.DataFrame()

        for layers in layer_iter:
            ids = layers.index \
                        .get_level_values("aln_id") \
                        .unique() \
                        .difference(aln_leftovers.index) \
                        .to_numpy()
            if len(ids) > 0:
                alignments = db0.query_alignments(aln_id=ids)
            else:
                alignments = pd.DataFrame()

            alignments = pd.concat([aln_leftovers, alignments])
            layers = pd.concat([layer_leftovers, layers])

            aln_end = alignments["read_id"] == alignments["read_id"].iloc[-1]
            aln_leftovers = alignments.loc[aln_end]
            alignments = alignments.loc[~aln_end]

            layer_end = layers.index.get_level_values("aln_id").isin(aln_leftovers.index)
            layer_leftovers = layers.loc[layer_end]
            layers = layers.loc[~layer_end]
            layer_alns = layers.index.get_level_values("aln_id")
            
            for ref_name,ref_alns in alignments.groupby("ref_name"):
                coords = self._alns_to_coords(ref_alns)
                cache = self._tables_to_tracks(coords, ref_alns, layers)
                for ret in cache.iter_reads_slice():
                    yield ret

        if len(aln_leftovers) > 0:
            for ref_name,ref_alns in aln_leftovers.groupby("ref_name"):
                coords = self._alns_to_coords(ref_alns)
                cache = self._tables_to_tracks(coords, ref_alns, layer_leftovers)
                for ret in cache.iter_reads_slice():
                    yield ret

    def _tables_to_tracks(self, coords, alignments, layers):
        tracks = dict()

        layer_alns = layers.index.get_level_values("aln_id")

        if self.prms.shared_refs_only or self.prms.min_coverage > 1:
            track_covs = alignments.loc[layer_alns, ["track_id"]] \
                                   .set_index(layers.index) \
                                   .reset_index("aln_id", drop=True) \
                                   .set_index("track_id", append=True) \
                                   .index.value_counts()

            mask = track_covs >= self.prms.min_coverage
            if not np.any(mask):
                mrefs = pd.Index([])
            elif self.prms.shared_refs_only:
                track_counts = pd.MultiIndex.from_tuples(track_covs[mask].index) \
                                   .get_level_values(0) \
                                   .value_counts()

                mrefs = track_counts.index[track_counts == len(self.alns)]
            else:
                mrefs = track_covs[mask].index.get_level_values(0).unique()

            
            layers = layers.loc[(mrefs,slice(None))]
            layer_alns = layers.index.get_level_values("aln_id")
            alignments = alignments.loc[layer_alns.unique()]

        aln_groups = alignments.groupby("track_id").groups
        for parent in self.alns:
            if parent.id in aln_groups:
                track_alns = alignments.loc[aln_groups[parent.id]]
                track_layers = layers.loc[layer_alns.isin(track_alns.index)]
            else:
                track_alns = alignments.iloc[:0] 
                track_layers = layers.iloc[:0]   

            track = AlnTrack(parent, coords, track_alns, track_layers)

            #if not track.empty:
            track.calc_layers(self.fn_layers)

            tracks[parent.name] = track

        tracks = Tracks(self, coords, tracks)

        if not tracks.all_empty:
            tracks.load_compare(alignments.index.to_numpy())
        tracks.calc_refstats()

        return tracks

    def _alns_to_coords(self, alns):
        ref_coord = RefCoord(
            alns["ref_name"].iloc[0],
            alns["ref_start"].min(),
            alns["ref_end"].max())
        return self.index.get_coord_space(
            ref_coord, self.conf.is_rna, load_kmers=True, kmer_trim=True)

    def iter_reads_old(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if read_filter is None:
            read_filter = self.prms.read_filter
        if max_reads is None:
            max_reads = self.prms.max_reads

        dbfile0,db0 = list(self.dbs.items())[0]

        aln_iter = db0.query_alignments(
            self._aln_track_ids,
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
                    yield (read_id, self.alns)

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

        layers = db.query_layers(self.db_layers, self._aln_track_ids, self.coords, ids)

        for track in self.alns:
            track_alns = alns[alns["track_id"] == track.id]

            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i].dropna(axis=1, how="all") #.set_index(self.layer_refs)

            if len(track_alns) > 0:
                name = track_alns["ref_name"].iloc[0]
                fwd = track_alns["fwd"].iloc[0]
                start = track_alns["ref_start"].min()
                end = track_alns["ref_end"].max()
                ref_coord = RefCoord(name, start, end, fwd)
                track_coords = self.index.get_coord_space(
                    ref_coord, self.conf.is_rna, load_kmers=True, kmer_trim=True)
            else:
                track_coords = None

            track.set_data(track_coords, track_alns, track_layers)

            if not track.empty:
                track.calc_layers(self.fn_layers)

        if not self.all_empty:
            self.load_compare(ids)

    
    def close(self):
        for filename, db in self.dbs.items():
            db.close()
        if self.prms.output_format != "db":
            self.io.close()
