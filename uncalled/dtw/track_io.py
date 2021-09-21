#!/usr/bin/env python3

import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections

from . import RefCoord
from .track2 import AlnTrack2
from ..index import load_index
from ..pore_model import PoreModel
from .. import config

class TrackIOParams(config.ParamGroup):
    _name = "track_io"
TrackIOParams._def_params(
    ("input", None, None, "Input track(s)"),
    ("output", None, None,  "Output track"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("index_prefix", None, str, "BWA index prefix"),
    ("load_mat", True, bool, "If true will load a matrix containing specified layers from all reads overlapping reference bounds"),
    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
#    ("layers", DEFAULT_LAYERS, list, "Layers to load"),
    #("mode", "r", str, "Read (r) or write (w) mode"),
    ignore_toml={"input", "output", "ref_bounds"}
)

class TrackIO:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("track_io", *args, **kwargs)

        self.dbs = dict()

        self.tracks = dict()
        self.track_dbs = dict()
        self.track_ids = dict()

        self.input_tracks = list()
        self.output_tracks = dict()
        self.prev_fast5 = dict()
        self.prev_read = dict()
        self.prev_aln = dict()

        self._load_dbs(self.prms.output, True)
        self._load_dbs(self.prms.input, False)
        self.input_tracks = tuple(self.input_tracks)

        self.index = load_index(self.prms.index_prefix)
        self.model = PoreModel(self.conf.pore_model) 
        self._set_ref_bounds(self.prms.ref_bounds)

    def _set_ref_bounds(self, ref_bounds):
        if ref_bounds is not None:
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
        #TODO sort (reindex?) DF my specified track_names to maintain order
        for i,row in db.query_track(track_names).iterrows():
            self._init_track(db, row.name)

            conf = config.Config(toml=row.config)
            t = AlnTrack2(db, row.id, row.name, row.desc, row.groups, conf)
            self.input_tracks.append(t)

            self.conf.load_config(conf)
            #if self.prms.index_prefix is None:
            #    self.prms.index_prefix = t.conf.track_io.index_prefix

    def _init_output_tracks(self, db, track_names):
        if track_names is None:
            name = os.path.splitext(os.path.basename(db.filename))[0]
            track_names = [name]

        for name in track_names:
            self._init_track(db, name)

            self.prev_fast5[name] = (None, None)
            self.prev_read[name] = None
            self.prev_aln[name] = -1

            track = AlnTrack2(db, None, name, name, "dtw", self.conf)
            self.output_tracks[name] = track
            db.init_track2(track)

    def init_alignment(self, read_id, fast5, group_name, layers, track_name=None):
        if track_name is None:
            if len(self.output_tracks) == 1:
                track_name, = self.output_tracks
            else:
                raise ValueError("Must specify track name when using multiple output tracks")

        track = self.output_tracks[track_name]
        db = track.db

        self.prev_aln[track_name] += 1
        aln_id = self.prev_aln[track_name]

        if fast5 == self.prev_fast5[track_name][0]:
            fast5_id = self.prev_fast5[track_name][1]
        else:
            fast5_id = db.init_fast5(fast5)
            self.prev_fast5[track_name] = (fast5, fast5_id)

        if self.prev_read[track_name] != read_id:
            db.init_read(read_id, fast5_id)
            self.prev_read[track_name] = read_id

        mref_start = layers.index.min()
        mref_end = layers.index.max()+1
        samp_start = layers["sample"].min()
        samp_end = layers["sample"].max()

        ref_bounds = self.index.mref_to_ref_bound(mref_start, mref_end, not self.conf.is_rna)

        track.alignments = pd.DataFrame({
                "id" : [aln_id],
                "track_id" : [track.id],
                "read_id" : [read_id],
                "ref_name" : [ref_bounds.ref_name],
                "ref_start" : [ref_bounds.start],
                "ref_end" : [ref_bounds.end],
                "fwd" :     [ref_bounds.fwd],
                "samp_start" : [samp_start],
                "samp_end" : [samp_end],
                "tags" : [""]}).set_index("id")

        track.db.init_alignment(track.alignments)

        track.layers = None
        track.add_layer_group(group_name, layers)

        track.coords = self.index.get_coord_space(track.aln_ref_coord(aln_id), self.conf.is_rna, kmer_shift=0)

        return track #aln_id

    def input_count(self):
        return len(self.input_tracks)

    def load_refs(self, ref_bounds=None, full_overlap=False):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if self.coords is None:
            raise ValueError("Must set ref bounds")

        for track in self.input_tracks:
            track.coords = self.coords
            alignments = track.db.query_alignments(track.id, coords=self.coords, full_overlap=full_overlap)

            if full_overlap:
                ids = self.alignments.index.to_numpy()
            else:
                ids = None

            layers = dict()
            for group in track.groups:
                layers[group] = track.db.query_layers("dtw", self.coords, ids)

            track.set_data(self.coords, alignments, layers)

        return self.input_tracks
    
    def close(self):
        for filename, db in self.dbs.items():
            db.close()
            
class TrackSQL:
    def __init__(self, sqlite_db):
        self.filename = sqlite_db
        new_file = not os.path.exists(sqlite_db)
        self.con = sqlite3.connect(sqlite_db)
        self.cur = self.con.cursor()

        if new_file:
            self.init_tables()

    def close(self):
        self.cur.execute("CREATE INDEX IF NOT EXISTS dtw_idx ON dtw (mref, aln_id);")
        self.cur.execute("CREATE INDEX IF NOT EXISTS bcaln_idx ON bcaln (mref, aln_id);")
        self.con.close()

    def init_tables(self):
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS track (
                id INTEGER PRIMARY KEY,
                name TEXT,
                desc TEXT,
                groups TEXT,
                config TEXT
            );""")
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

    def init_track2(self, track):
        self.cur.execute(
            "INSERT INTO track (name,desc,config,groups) VALUES (?,?,?,?)",
            (track.name, track.desc, track.conf.to_toml(), "dtw")
        )
        track.id = self.cur.lastrowid
        #self.con.commit()
        return track.id

    def init_track(self, track):
        self.cur.execute(
            "INSERT INTO track (name,desc,config,groups) VALUES (?,?,?,?)",
            (track.prms.name, track.prms.name, track.conf.to_toml(), "dtw")
            #(track.name, track.desc, track.config.to_toml(), "dtw")
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
        #self.con.commit()

    def write_layers(self, table, df):
        df.to_sql(
            table, self.con, 
            if_exists="append", 
            method="multi",# chunksize=5000,
            index=True, index_label=["mref","aln_id"])
        #self.con.commit()

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
        

    def query_alignments(self, track_id=None, read_id=None, coords=None, full_overlap=False):
        select = "SELECT * FROM alignment"
        wheres = list()
        params = list()

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

        query = self._join_query(select, wheres)
        
        return pd.read_sql_query(query, self.con, params=params).set_index("id")
        
    def _join_query(self, select, wheres):
        if len(wheres) == 0:
            return select
        else:
            return select + " WHERE " + " AND ".join(wheres)

    def query_layers(self, table, coords=None, aln_id=None, order="mref"):
        select = "SELECT mref, aln_id, start, length, current FROM dtw"

        wheres = list()
        params = list()

        if coords is not None:
            wheres.append("mref >= ? AND mref <= ?")
            params += [str(coords.mrefs[True].min()), str(coords.mrefs[True].max())]

        self._add_where(wheres, params, "aln_id", aln_id)

        #if len(wheres) == 0:
        #else:
        #    query = select + " WHERE " + " AND ".join(wheres)
        query = self._join_query(select, wheres)

        
        return pd.read_sql_query(query, self.con, params=params).set_index(["mref","aln_id"])
