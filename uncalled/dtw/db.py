"""Subcommands to alter DTW track databases

subcommand options:
ls        List all tracks in a database
delete    Delete a track from a database
edit     Rename and/or set the description of a track"""
#copy      Copy a track from one database to another"""

import sqlite3
import os
import collections
import numpy as np
import pandas as pd

from ..argparse import Opt

DB_OPT = Opt("db_file", help="Track database file")

LS_OPTS = (DB_OPT,)
_LS_QUERY = "SELECT name,desc,COUNT(alignment.id) FROM track " \
            "JOIN alignment ON track.id == track_id GROUP BY name"
def ls(conf, con=None):
    if con is None:
        con = sqlite3.connect(conf.db_file)
    cur = con.cursor()
    print("\t".join(["Name", "Description", "Alignments"]))
    
    for row in cur.execute(_LS_QUERY).fetchall():
        print("\t".join(map(str, row)))
    con.commit()

def _verify_track(cur, track_name):
    existing = cur.execute("SELECT name FROM track WHERE name == ?", (track_name,))
    if len(existing.fetchall()) == 0:
        sys.stderr.write("Track does not exist: \"%s\"\n" % conf.track_name)
        sys.exit(1)

DELETE_OPTS = (
    DB_OPT,
    Opt("track_name", help="Name of the track to delete"),
)
def delete(conf, con=None):
    if con is None:
        con = sqlite3.connect(conf.db_file)
    cur = con.cursor()
    _verify_track(cur, conf.track_name)

    cur.execute("PRAGMA foreign_keys = ON")
    cur.execute("DELETE FROM track WHERE name == ?", (conf.track_name,))
    con.commit()
    print("Deleted track \"%s\"" % conf.track_name)

EDIT_OPTS = (
    Opt("db_file", help="Track database file"),
    Opt("track_name", help="Current track name"),
    Opt(("-N", "--new-name"), default=None, help="New track name"),
    Opt(("-D", "--description"), default=None, help="New track description"),
)
def edit(conf, con=None):
    if con is None: con = sqlite3.connect(conf.db_file)
    cur = con.cursor()
    _verify_track(cur, conf.track_name)

    updates = []
    params = []
    if conf.new_name:
        updates.append("name = ?")
        params.append(conf.new_name)
    if conf.description:
        updates.append("desc = ?")
        params.append(conf.description)
    if len(updates) == 0:
        sys.stderr.write("No changes made. Must specify new name (-N) or description (-D)\n")
        return
    params.append(conf.track_name)
        
    query = "UPDATE track SET " + ", ".join(updates) + " WHERE name == ?"
    cur.execute(query, params)
    con.close()

SUBCMDS = [
    (ls, LS_OPTS), 
    (delete, DELETE_OPTS), 
    (edit, EDIT_OPTS), 
]


            
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
                FOREIGN KEY (track_id) REFERENCES track (id) ON DELETE CASCADE,
                FOREIGN KEY (read_id) REFERENCES read (id)
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS dtw (
                mref INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                current REAL,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS bcaln (
                mref INTEGER,
                aln_id INTEGER,
                sample INTEGER,
                bp INTEGER,
                error TEXT,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
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

    def get_fast5_index(self, track_id=None):
        query = "SELECT read.id AS read_id, filename FROM read " +\
                "JOIN fast5 ON fast5.id = fast5_id"

        wheres = list()
        params = list()

        if track_id is not None:
            query += " JOIN alignment ON read.id = alignment.read_id"
            self._add_where(wheres, params, "track_id", track_id)
        
        query = self._join_query(query, wheres)
        return pd.read_sql_query(query, self.con, params=params)

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
            if coords.stranded:
                wheres.append("(mref >= ? AND mref <= ?)")
                params += [str(coords.mrefs.min()), str(coords.mrefs.max())]
            else:
                wheres.append("((mref >= ? AND mref <= ?) OR (mref >= ? AND mref <= ?))")
                params += [str(coords.mrefs[True].min()), str(coords.mrefs[True].max())]
                params += [str(coords.mrefs[False].min()), str(coords.mrefs[False].max())]

        self._add_where(wheres, params, "aln_id", aln_id)

        query = self._join_query(select, wheres, order)

        return pd.read_sql_query(query, self.con, index_col=index, params=params, chunksize=chunksize)

        #return pd.concat({table : df}, names=["group", "layer"], axis=1)
