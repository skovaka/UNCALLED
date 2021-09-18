#!/usr/bin/env python3

import sys, os
import sqlite3
import numpy as np
import pandas as pd

from .. import config

#class TrackIOParams(config.ParamGroup):
#    _name = "track_io"
#TrackIOParams.def_params(
#    ("input", None, None, "Input track(s)"),
#    ("output", str, None, "Output track"),
#    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
#    ("index_prefix", None, str, "BWA index prefix"),
#    ("load_mat", True, bool, "If true will load a matrix containing specified layers from all reads overlapping reference bounds"),
#    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
#    ("layers", DEFAULT_LAYERS, list, "Layers to load"),
#    #("mode", "r", str, "Read (r) or write (w) mode"),
#    ignore_toml={"mode", "overwrite"}
#)
#
#class TrackIO:
#    def __init__(self, *args, **kwargs):
#        self.conf, self.prms = config._init_group("track_io", *args, **kwargs)
#        
#        self.dbs = list()
#
#        if self.prms.input is None:
#            self.in_dbs = None
#        else:
#            if not isinstance(self.prms.input, list):
#                self.prms.input = [self.prms.input]
#            self.in_dbs = list()
#            for db in self.prms.input:
#                track_names = self.open_db(db)
#                self.in_dbs.append()
#
#    def open_db(self, db_str):
#        spl = db_str.split(":")
#        if len(spl) == 2:
#            filename = spl[0]
#            track_names = spl[1].split(",")
#        elif len(spl) == 1:
#            filename = db_str
#            track_names = None
#        else:
#            raise ValueError("Incorrect database specifier format: " + db_str)
#
#        db = TrackSQL(filename)


class TrackSQL:
    def __init__(self, sqlite_db):
        self.con = sqlite3.connect(sqlite_db)
        self.cur = self.con.cursor()
        #self.prms = track_io_prms

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
                config TEXT,
                groups TEXT
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

        #self.cur.execute(
        #    "INSERT INTO alignment (track_id, read_id, ref_name, ref_start, ref_end, fwd) VALUES (?,?,?,?,?,?)",
        #    (aln.track_id, aln.read_id, aln.ref_name, aln.ref_start, aln.ref_end, aln.is_fwd)
        #)
        #aln.id = self.cur.lastrowid
        #self.con.commit()

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

    def query_track(self, name):
        return self.cur.execute(
            "SELECT id, desc, config, groups FROM track WHERE name == ?", 
            (name,)).fetchone()
        
    def query_read(self, read_id, coords=None):
        query = """
            SELECT mref, aln_id, start, length, current FROM dtw
            JOIN alignment ON id = aln_id
            WHERE read_id = ?"""
        params = [read_id]

        if coords is not None:
            query = select + " AND (mref >= ? AND mref <= ?)"
            params += [int(coords.mrefs[True].min()), str(coords.mrefs[True].max())]

        dtw = pd.read_sql_query(query, self.con, params=params).set_index("mref")
        aln_id = str(dtw["aln_id"].iloc[0])

        aln = pd.read_sql_query(
            "SELECT * FROM alignment WHERE id = ?",
            self.con, params=(aln_id,)
        ).set_index("id")

        return aln, dtw

    def query_alns(self, track_id, coords, full_overlap):

        select = "SELECT mref, aln_id, start, length, current FROM dtw"
        where = " WHERE (mref >= ? AND mref <= ?)"
        params = [int(coords.mrefs[True].min()), str(coords.mrefs[True].max())]

        if full_overlap:
            select = select + " JOIN alignment ON id = aln_id"
            where = where + " AND (ref_start < ? AND ref_end > ?)"
            params += [int(coords.refs.min()), int(coords.refs.max())]

        query = select + where

        dtw = pd.read_sql_query(query, self.con, params=params)

        ids = dtw["aln_id"][~dtw["aln_id"].duplicated()].to_numpy(dtype=str)

        alns = pd.read_sql_query(
            "SELECT * FROM alignment WHERE id IN (%s)" % ",".join(["?"]*len(ids)),
            self.con, params=ids
        ).sort_values("ref_start")

        return alns, dtw
