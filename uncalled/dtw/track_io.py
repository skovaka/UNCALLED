#!/usr/bin/env python3

import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections
import time

from .track import AlnTrack, DEFAULT_LAYERS
from ..index import load_index, RefCoord, str_to_coord
from ..pore_model import PoreModel
from ..fast5 import Fast5Reader, parse_read_ids
from .. import config

class TrackIOParams(config.ParamGroup):
    _name = "track_io"
TrackIOParams._def_params(
    ("input", None, None, "Input track(s)"),
    ("output", None, None,  "Output track"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("index_prefix", None, str, "BWA index prefix"),
    ("overwrite", False, bool, "Overwrite existing databases"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("aln_chunksize", 4000, int, "Number of alignments to query for iteration"),
    ("ref_chunksize", 10000, int, "Number of reference coordinates to query for iteration"),
    ("layers", ["current", "dwell", "model_diff"], list, "Layers to load"),
    #("mode", "r", str, "Read (r) or write (w) mode"),
    ignore_toml={"input", "output", "ref_bounds", "layers"}
)

class TrackIO:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("track_io", *args, **kwargs)

        self.dbs = dict()

        self.track_dbs = dict()

        self.input_tracks = list()
        self.output_tracks = dict()
        self.prev_fast5 = dict()
        self.prev_read = dict()

        self._load_dbs(self.prms.output, True)
        self._load_dbs(self.prms.input, False)

        self.input_track_ids = [t.id for t in self.input_tracks]
        self.output_track_ids = [t.id for _,t in self.output_tracks.items()]

        if self.prms.index_prefix is not None:
            self.index = load_index(self.prms.index_prefix)
            self._set_ref_bounds(self.prms.ref_bounds)

        self.model = PoreModel(self.conf.pore_model) 
        self.fast5s = None

    def _set_ref_bounds(self, ref_bounds):
        if ref_bounds is not None:
            if isinstance(ref_bounds, str):
                ref_bounds = str_to_coord(ref_bounds)
            elif isinstance(ref_bounds, tuple):
                ref_bounds = RefCoord(*ref_bounds)
            self.coords = self.index.get_coord_space(ref_bounds, self.conf.is_rna)
        else:
            self.coords = None

    def _load_dbs(self, dbs, out):
        if dbs is None:
            return

        if isinstance(dbs, str):
            dbs = [dbs]

        for db_str in dbs:
            db_file, track_names = self._db_track_split(db_str)

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
        #TODO sort (reindex?) DF by specified track_names to maintain order
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
            self.input_tracks.append(t)

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
        return len(self.input_tracks)


    LAYER_FNS = {
        "dwell" : (lambda self,track: 
            1000 * track.layers["dtw","length"] / self.conf.read_buffer.sample_rate),
        "model_diff" : (lambda self,track: 
            track.layers["dtw","current"] - self.model[track.kmers])
    }

    def compute_layers(self, track, layer_names):
        for name in layer_names:
            if not name in track.layers.columns.get_level_values(1):
                track.layers["dtw",name] = self.LAYER_FNS[name](self,track)

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
        
        for track in self.input_tracks:
            track_alns = alignments[alignments["track_id"] == track.id].copy()
            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i].copy()

            track.set_data(self.coords, track_alns, track_layers)
            self.compute_layers(track, self.prms.layers)
            
            if load_mat:
                track.load_mat()

        return self.input_tracks

    def iter_refs(self, ref_bounds=None):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)

        dbfile0,db0 = list(self.dbs.items())[0]
        layer_iter = db0.query_layers(
            self.input_track_ids, 
            coords=self.coords, 
            order=["mref"],
            chunksize=self.prms.ref_chunksize)

        end_layers = None
        seq_coords = None
        for chunk in layer_iter:
            if end_layers is not None:
                chunk = pd.concat([end_layers, chunk])

            chunk_mrefs = chunk.index.get_level_values("mref").unique()

            if seq_coords is not None:
                coords = seq_coords.mref_intersect(chunk_mrefs[:-1])
                
            if seq_coords is None or coords is None:
                chunk_refs = self.index.mrefs_to_ref_coord(chunk_mrefs[0], chunk_mrefs[-1], not self.conf.is_rna)
                seq_refs = RefCoord(chunk_refs.name, 0, chunk_refs.ref_len, chunk_refs.fwd)
                seq_coords = self.index.get_coord_space(seq_refs, self.conf.is_rna, kmer_shift=0, load_kmers=False)
                coords = seq_coords.mref_intersect(chunk_mrefs[:-1])

            coords.kmers = self.index.mrefs_to_kmers(coords.mrefs, self.conf.is_rna)

            leftover = chunk_mrefs.difference(coords.mrefs)
            end_layers = chunk.loc[leftover]
            chunk = chunk.drop(index=leftover)

            layers = pd.concat({"dtw":chunk}, names=["group", "layer"], axis=1)

            aln_ids = chunk.index.unique("aln_id").to_numpy()
            alns = db0.query_alignments(self.input_track_ids, aln_id=aln_ids)

            for track in self.input_tracks:
                track_alns = alns[alns["track_id"]==track.id].copy()
                track_layers = layers[layers.index.isin(track_alns.index, 1)].copy()
                track.set_data(coords, track_alns, track_layers)
                self.compute_layers(track, self.prms.layers)

            yield (coords, self.input_tracks)

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

        return self.input_tracks

    def iter_reads(self, ref_bounds=None, full_overlap=False, max_reads=None):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        
        dbfile0,db0 = list(self.dbs.items())[0]

        aln_iter = db0.query_alignments(
            self.input_track_ids,
            coords=self.coords, 
            full_overlap=full_overlap, 
            order=["read_id"],
            chunksize=self.prms.aln_chunksize)

        n = 0

        end_alns = None
        for chunk in aln_iter:
            if end_alns is not None:
                chunk = pd.concat([end_alns, chunk])

            end = chunk["read_id"] == chunk["read_id"].iloc[-1]
            end_alns = chunk[end]
            chunk = chunk[~end]

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
                    n += 1
                    if max_reads is not None and n >= max_reads:
                        return (read_id, self.input_tracks)
                    yield (read_id, self.input_tracks)

    def _fill_tracks(self, db, alns):
        ids = list(alns.index)

        layers = dict()
        #TODO parse which track layers to import
        for group in ["dtw"]: #self.groups:
            layers[group] = db.query_layers(self.input_track_ids, self.coords, ids, index=["aln_id","mref"])
        layers = pd.concat(layers, names=["group", "layer"], axis=1)

        for track in self.input_tracks:
            track_alns = alns[alns["track_id"] == track.id]

            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i]

            name = track_alns["ref_name"].iloc[0]
            fwd = track_alns["fwd"].iloc[0]
            start = track_alns["ref_start"].min()
            end = track_alns["ref_end"].max()
            ref_coord = RefCoord(name, start, end, fwd)
            track_coords = self.index.get_coord_space(
                ref_coord, self.conf.is_rna, kmer_shift=0, load_kmers=True)

            track.set_data(track_coords, track_alns, track_layers)
            self.compute_layers(track, self.prms.layers)
    
    def close(self):
        for filename, db in self.dbs.items():
            db.close()
            
class TrackSQL:
    def __init__(self, sqlite_db):
        self.filename = sqlite_db
        new_file = not os.path.exists(sqlite_db)
        self.con = sqlite3.connect(sqlite_db)
        self.cur = self.con.cursor()
        self.open = True

        #if new_file:
        self.init_tables()

        self.prev_aln_id = self.cur.execute("SELECT MAX(id) FROM alignment").fetchone()[0]
        if self.prev_aln_id is None:
            self.prev_aln_id = -1

    def close(self):
        if self.open:
            self.cur.execute("CREATE INDEX IF NOT EXISTS dtw_idx ON dtw (mref, aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS bcaln_idx ON bcaln (mref, aln_id);")
            self.con.close()
            self.open = False

    def init_tables(self):
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS track (
                id INTEGER PRIMARY KEY,
                name TEXT,
                desc TEXT,
                groups TEXT,
                config TEXT
            );""")
        self.cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS track_idx ON track (name);")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS read (
                id TEXT PRIMARY KEY,
                fast5_id INTEGER,
                FOREIGN KEY (fast5_id) REFERENCES fast5 (id)
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS fast5 (
                id INTEGER PRIMARY KEY,
                filename TEXT
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS alignment (
                id INTEGER PRIMARY KEY,
                track_id INTEGER,
                read_id TEXT,
                ref_name TEXT,
                ref_start INTEGER,
                ref_end INTEGER,
                fwd INTEGER,
                samp_start INTEGER,
                samp_end INTEGER,
                tags TEXT,
                FOREIGN KEY (track_id) REFERENCES track (id),
                FOREIGN KEY (read_id) REFERENCES read (id)
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS dtw (
                mref INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                current REAL,
                FOREIGN KEY (aln_id) REFERENCES alignment (id)
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS bcaln (
                mref INTEGER,
                aln_id INTEGER,
                sample INTEGER,
                bp INTEGER,
                error TEXT,
                FOREIGN KEY (aln_id) REFERENCES alignment (id)
            );""")
        #self.con.commit()

    def init_write(self):
        for table in ["dtw", "bcaln"]:
            self.cur.execute("DROP INDEX IF EXISTS %s_idx" % table)

    def init_track(self, track):
        self.cur.execute(
            "INSERT INTO track (name,desc,config,groups) VALUES (?,?,?,?)",
            (track.name, track.desc, track.conf.to_toml(), "dtw")
        )
        track.id = self.cur.lastrowid
        #self.con.commit()
        return track.id

    def init_alignment(self, aln_df):
        aln_df.to_sql(
            "alignment", self.con, 
            if_exists="append", 
            method="multi",# chunksize=5000,
            index=True, index_label="id")

    def init_fast5(self, filename):

        row = self.cur.execute("SELECT id FROM fast5 WHERE filename = ?", (filename,)).fetchone()
        if row is not None:
            return row[0]
            
        self.cur.execute("INSERT INTO fast5 (filename) VALUES (?)", (filename,))
        fast5_id = self.cur.lastrowid
        #self.con.commit()

        return fast5_id

    def init_read(self, read_id, fast5_id):
        self.cur.execute("INSERT OR IGNORE INTO read VALUES (?,?)", (read_id, fast5_id))

    def next_aln_id(self):
        self.prev_aln_id += 1
        return self.prev_aln_id

    def write_layers(self, df):
        for group in df.columns.levels[0]:
            df[group].to_sql(
                group, self.con, 
                if_exists="append", 
                method="multi", chunksize=50000,
                index=True, index_label=["mref","aln_id"])

    def get_fast5_index(self):
        return pd.read_sql_query("""
            SELECT read.id AS read_id, filename FROM read
            JOIN fast5 ON fast5.id = fast5_id""", self.con)

    def get_read_ids(self):
        return set(pd.read_sql_query("SELECT id FROM read", self.con)["id"])

    def _add_where(self, wheres, params, name, val):
        if val is None: return

        if not isinstance(val, str) and isinstance(val, (collections.abc.Sequence, np.ndarray)):
            wheres.append("%s in (%s)" % (name, ",".join(["?"]*len(val))))
            params += map(str, val)
        else:
            wheres.append("%s = ?" % name)
            params.append(str(val))

    def query_track(self, name=None):
        select = "SELECT * FROM track"
        wheres = list()
        params = list()
        self._add_where(wheres,params,"name",name)
        query = self._join_query(select, wheres)
        return pd.read_sql_query(query, self.con, params=params)
        #return self.cur.execute(query,params).fetchone()

    def query_alignments(self, track_id=None, read_id=None, aln_id=None, coords=None, full_overlap=False, order=None, chunksize=None):
        select = "SELECT * FROM alignment"
        wheres = list()
        params = list()

        self._add_where(wheres, params, "id", aln_id)
        self._add_where(wheres, params, "track_id", track_id)
        self._add_where(wheres, params, "read_id", read_id)

        if coords is not None:
            ref_start = int(coords.refs.min())
            ref_end = int(coords.refs.max())+1
            wheres.append("ref_name = ? AND ref_start < ? AND ref_end > ?")
            params.append(coords.ref_name)
            if full_overlap:
                params += [ref_start, ref_end]
            else:
                params += [ref_end, ref_start]

        query = self._join_query(select, wheres, order)
        
        return pd.read_sql_query(query, self.con, index_col="id", params=params, chunksize=chunksize)
        
    def _join_query(self, select, wheres, order=None):
        query = select
        if len(wheres) > 0:
            query += " WHERE " + " AND ".join(wheres)
        if order is not None:
            query += " ORDER BY " + ", ".join(order)
        return query

    def query_layers(self, track_id=None, coords=None, aln_id=None, order=["mref"], index=["mref","aln_id"], chunksize=None):
        select = "SELECT mref, aln_id, start, length, current FROM dtw"

        wheres = list()
        params = list()

        if track_id is not None:
            select += " JOIN alignment ON id = aln_id"
            self._add_where(wheres, params, "track_id", track_id)

        if coords is not None:
            wheres.append("((mref >= ? AND mref <= ?) OR (mref >= ? AND mref <= ?))")
            params += [str(coords.mrefs[True].min()), str(coords.mrefs[True].max())]
            params += [str(coords.mrefs[False].min()), str(coords.mrefs[False].max())]

        self._add_where(wheres, params, "aln_id", aln_id)

        query = self._join_query(select, wheres, order)

        return pd.read_sql_query(query, self.con, index_col=index, params=params, chunksize=chunksize)

        #return pd.concat({table : df}, names=["group", "layer"], axis=1)
