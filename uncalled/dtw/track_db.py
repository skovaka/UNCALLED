#!/usr/bin/env python3

import sys, os
import numpy as np
import sqlite3

#class TrackIOParams(ParamGroup):
#    _name = "track_io"
#TrackIOParams.def_params(
#    ("path", None, None, "Path to directory where alignments are stored"),
#)

class TrackSQL:
    def __init__(self, sqlite_db):
        self.con = sqlite3.connect(self.db_filename)
        #self.prms = track_io_prms

    def init_tables(self):
        cur = self.con.cursor()
        cur.execute("""
            CREATE TABLE IF NOT EXISTS track (
                id INTEGER PRIMARY KEY,
                name TEXT,
                desc TEXT,
                config TEXT,
                groups TEXT
            );
            CREATE TABLE IF NOT EXISTS read (
                id TEXT PRIMARY KEY,
                name TEXT,
                desc TEXT,
                config TEXT
            );
            CREATE TABLE IF NOT EXISTS fast5 (
                id INTEGER PRIMARY KEY,
                filename TEXT
            );
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
            );
            CREATE TABLE IF NOT EXISTS dtw (
                mref INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                current REAL,
                PRIMARY KEY (mref, aln_id),
                FOREIGN KEY (aln_id) REFERENCES alignment (id)
            );
            CREATE TABLE IF NOT EXISTS bcaln (
                mref INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                error TEXT,
                PRIMARY KEY (mref, aln_id),
                FOREIGN KEY (aln_id) REFERENCES alignment (id)
            );
            """

    def write_track(self, track):
        cur = self.con.cursor()
        cur.execute(
            "INSERT INTO track (name,desc,config,groups) VALUES (?,?,?,?)",
            (track.name, track.desc, track.config.to_toml(), "dtw")
        )
        self.con.commit()

    def write_aln(self, aln):
        cur = self.con.cursor()
        samp_st, samp_end = get_samp_bounds(self)
        cur.execute(
            "INSERT INTO alignment (track_id, read_id, ref_name, ref_start, ref_end, fwd, samp_start, samp_end) VALUES (?,?,?,?,?,?,?,?)",
            (aln.track_id, aln.read_id, aln.ref_name, aln.ref_start, aln.ref_end, aln.is_fwd, samp_start, samp_end)
        )
        aln.id = cur.lastrowid
        self.write_layers("dtw", aln.dtw)
        self.con.commit()

    def write_layers(self, table, df):
        df.to_sql(table, self.con, if_exists="append", index=True, index_label="mref")

    def query(self, track_id, mref_coords, full_overlap):

        select = "SELECT mref, aln_id, start, length, current FROM dtw"
        where = " WHERE (mref >= ? AND mref <= ?)"
        params = [int(mref_coords[True].min()), str(mref_coords[True].max())]

        if full_overlap:
            select = select + " JOIN alignment ON id = aln_id"
            where = where + " AND (ref_start < ? AND ref_end > ?)"
            params += [int(mref_coords.index[0]), int(mref_coords.index[-1])]

        query = select + where

        dtw = pd.read_sql_query(query, self.con, params=params)

        ids = dtw["aln_id"][~dtw["aln_id"].duplicated()].to_numpy(dtype=str)

        alns = pd.read_sql_query(
            "SELECT * FROM alignment WHERE id IN (%s)" % ",".join(["?"]*len(ids),
            self.con, params=ids
        ).sort_values("ref_start")

        return alns, dtw
