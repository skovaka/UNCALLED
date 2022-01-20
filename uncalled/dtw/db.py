"""Edit, merge, and ls alignment databases

subcommand options:
ls       List all tracks in a database
delete   Delete a track from a database
merge    Merge databases into a single file
edit     Rename, change fast5 paths, or set description"""

import sqlite3
import os
import collections
import numpy as np
import pandas as pd
import sys
import time

from ..nt import kmer_to_str, kmer_rev, kmer_comp
from ..argparse import Opt, comma_split
from ..fast5 import parse_fast5_paths

class TrackIO:
    def __init__(self):
        if not hasattr(self, "prev_aln_id"):
            self.prev_aln_id = 0

    def next_aln_id(self):
        self.prev_aln_id += 1
        return self.prev_aln_id


class Eventalign(TrackIO):
    def __init__(self, filename=None, mode="r"):
        self.mode = "r"
        if filename is None:
            self.out = sys.stdout
        else:
            self.out = open(filename, mode)

        self._header = True

        TrackIO.__init__(self)

    def write_events(self, track, events):
        contig = track.coords.ref_name

        events = events.set_index("mref")
        kmers = track.coords.kmers[events.index]
        model_kmers = kmer_to_str(kmer_rev(kmers))

        if track.coords.fwd:
            ref_kmers = model_kmers
        else:
            ref_kmers = kmer_to_str(kmer_comp(kmers))

        #read_id = track.alignments.iloc[0]["read_id"]
        #sys.stderr.write(f"{self.prev_aln_id}\t{read_id}\n")

        model = track.model

        #https://github.com/jts/nanopolish/issues/655
        std_level = (events["current"] - model.model_mean) / model.model_stdv

        eventalign = pd.DataFrame(
            data = {
                "contig" : track.coords.ref_name,
                "position" : track.coords.mref_to_ref(events.index)-2,
                "reference_kmer" : ref_kmers,
                "read_index" : self.prev_aln_id,
                #"read_name" : read_id,
                "strand" : "t",
                "event_index" : pd.RangeIndex(0,len(events))[::-1]+1,
                "event_level_mean" : events["current"],
                "event_stdv" : events["current_stdv"],
                "event_length" : events["length"] / track.conf.read_buffer.sample_rate,
                "model_kmer" : model_kmers,
                "model_mean" : model.means[kmers],
                "model_stdv" : model.stdvs[kmers],
                "standardized_level" : std_level,
                "start_idx" : events["start"],
                "end_idx" : events["start"] + events["length"],
            }, index = events.index).sort_values("position")

        eventalign.to_csv(
            self.out, sep="\t", 
            header=self._header,
            float_format="%.5f",
            index=False)

        self._header = False

    def init_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        self.out.close()


class TrackSQL(TrackIO):
    def __init__(self, sqlite_db):
        self.filename = sqlite_db
        new_file = not os.path.exists(sqlite_db)
        self.con = sqlite3.connect(sqlite_db)
        self.cur = self.con.cursor()
        self.cur.execute("PRAGMA foreign_keys = ON")
        self.con.commit()
        self.open = True

        #if new_file:
        self.init_tables()

        self.prev_aln_id = self.cur.execute("SELECT MAX(id) FROM alignment").fetchone()[0]
        if self.prev_aln_id is None:
            self.prev_aln_id = 0

        TrackIO.__init__(self)

    def close(self):
        if self.open:
            self.cur.execute("CREATE INDEX IF NOT EXISTS dtw_idx ON dtw (mref, aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS bcaln_idx ON bcaln (mref, aln_id);")
            #self.cur.execute("CREATE INDEX IF NOT EXISTS cmp_idx ON cmp (mref, aln_id, aln_b, group_b);")
            self.con.close()
            self.open = False

    def init_tables(self):
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS track (
                id INTEGER PRIMARY KEY,
                name TEXT,
                desc TEXT,
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
                id INTEGER PRIMARY KEY ON CONFLICT REPLACE,
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
                stdv REAL,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS bcaln (
                mref INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                bp INTEGER,
                error TEXT,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS cmp (
                mref INTEGER,
                aln_a INTEGER,
                aln_b INTEGER,
                group_b TEXT DEFAULT "dtw",
                mean_ref_dist REAL, 
                jaccard REAL, 
                FOREIGN KEY (aln_a) REFERENCES alignment (id) ON DELETE CASCADE,
                FOREIGN KEY (aln_b) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        #self.con.commit()

    def init_write(self):
        for table in ["dtw", "bcaln", "cmp"]:
            self.cur.execute("DROP INDEX IF EXISTS %s_idx" % table)

    def init_track(self, track):
        self.cur.execute(
            "INSERT INTO track (name,desc,config) VALUES (?,?,?)",
            (track.name, track.desc, 
             track.conf.to_toml())
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

    def write_layers(self, df, index=["mref","aln_id"]):
        for group in df.columns.levels[0]:
            df[group].to_sql(
                group, self.con, 
                if_exists="append", 
                method="multi", chunksize=50000,
                index=True, index_label=index)

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
            if len(val) == 0: return
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

    def query_compare(self, layers, track_id=None, coords=None, aln_id=None):
        dtw = False
        bcaln = False
        fields = {"mref", "aln_a", "aln_b", "group_b"}
        for group,layer in layers:
            if group == "cmp": dtw = True
            if group == "bc_cmp": bcaln = True
            fields.add(layer)

        fields = ", ".join(fields)
        select = f"SELECT {fields} FROM cmp"

        wheres = list()
        params = list()

        if track_id is not None:
            select += " JOIN alignment ON id = cmp.aln_a"
            self._add_where(wheres, params, "track_id", track_id)

        if dtw and not bcaln:
            self._add_where(wheres, params, "group_b", "dtw")
        elif bcaln and not dtw:
            self._add_where(wheres, params, "group_b", "bcaln")

        self._add_where(wheres, params, "aln_a", aln_id)

        if coords is not None:
            if coords.stranded:
                wheres.append("(mref >= ? AND mref <= ?)")
                params += [str(coords.mrefs.min()), str(coords.mrefs.max())]
            else:
                wheres.append("((mref >= ? AND mref <= ?) OR (mref >= ? AND mref <= ?))")

                params += [str(coords.mrefs[True].min()), str(coords.mrefs[True].max())]
                params += [str(coords.mrefs[False].min()), str(coords.mrefs[False].max())]

        query = self._join_query(select, wheres)

        return pd.read_sql_query(
            query, self.con, 
            index_col=["mref", "aln_a", "aln_b", "group_b"], 
            params=params)


    def query_layers(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, order=["mref"], chunksize=None, full_overlap=False):

        t0 = time.time()

        group_layers = collections.defaultdict(list)
        renames = dict()
        fields = list()
        tables = list()
        for group,layer in layers:
            if group == "bc_cmp":
                group = "cmp"
            field = group + "." + layer
            name = group + "_" + layer
            group_layers[group].append(name)
            renames[name] = layer
            fields.append(field + " AS " + name)
            if group not in tables:
                tables.append(group)

        fields += ["%s.%s AS idx_%s" % (tables[0],idx,idx) for idx in ("mref", "aln_id")]

        select = "SELECT " + ", ".join(fields) + " FROM " + tables[0]
        for table in tables[1:]:
            select += " LEFT JOIN %s ON %s.aln_id == idx_aln_id AND %s.mref == idx_mref" % ((table,)*3)

        wheres = list()
        params = list()

        if track_id is not None or read_id is not None:
            select += " JOIN alignment ON id = idx_aln_id"
            if track_id is not None:
                self._add_where(wheres, params, "track_id", track_id)
            if read_id is not None:
                self._add_where(wheres, params, "read_id", read_id)

        if "cmp" in group_layers:
            self._add_where(wheres, params, "group_b", "bcaln")

        if coords is not None:
            if coords.stranded:
                wheres.append("(idx_mref >= ? AND idx_mref <= ?)")
                params += [str(coords.mrefs.min()), str(coords.mrefs.max())]
            else:
                wheres.append("((idx_mref >= ? AND idx_mref <= ?) OR (idx_mref >= ? AND idx_mref <= ?))")

                params += [str(coords.mrefs[True].min()), str(coords.mrefs[True].max())]
                params += [str(coords.mrefs[False].min()), str(coords.mrefs[False].max())]
            if full_overlap:
                wheres.append("ref_name = ? AND ref_start < ? AND ref_end > ?")
                params.append(coords.ref_name)
                params += [coords.refs.min(), coords.refs.max()]

        self._add_where(wheres, params, "idx_aln_id", aln_id)

        query = self._join_query(select, wheres, ["idx_"+o if o in {"aln_id","mref"} else o for o in order])

        t0 = time.time()

        ret = pd.read_sql_query(
            query, self.con, 
            index_col=["idx_mref", "idx_aln_id"], 
            params=params, chunksize=chunksize)

        t0 = time.time()

        def make_groups(df):
            grouped = dict()
            for group, layers in group_layers.items():
                gdf = df[layers]
                grouped[group] = df[layers].rename(columns=renames)
            df = pd.concat(grouped, names=("group", "layer"), axis=1)
            df.index.names = ("mref", "aln_id")
            return df
                
        if chunksize is None:
            return make_groups(ret)

        t0 = time.time()

        return (make_groups(df) for df in ret)
    
    def _verify_track(self, track_name):
        ids = self.cur.execute("SELECT id FROM track WHERE name == ?", (track_name,)).fetchall()
        if len(ids) == 0:
            raise ValueError(f"Track does not exist: \"{track_name}\"\n")
        return ids[0][0]

DB_OPT = Opt("db_file", help="Track database file")

LS_OPTS = (DB_OPT,)
_LS_QUERY = "SELECT name,desc,COUNT(alignment.id) FROM track " \
            "JOIN alignment ON track.id == track_id GROUP BY name"
def ls(conf, db=None):
    if db is None:
        db = TrackSQL(conf.db_file)
    print("\t".join(["Name", "Description", "Alignments"]))

    for row in db.cur.execute(_LS_QUERY).fetchall():
        print("\t".join(map(str, row)))

    db.con.commit()

DELETE_OPTS = (
    DB_OPT,
    Opt("track_name", help="Name of the track to delete"),
)
def delete(track_name=None, db=None, conf=None):
    if db is None:
        db = TrackSQL(conf.db_file)

    if track_name is None:
        track_name = conf.track_name

    db._verify_track(track_name)

    db.cur.execute("DELETE FROM track WHERE name == ?", (track_name,))
    db.con.commit()
    print("Deleted track \"%s\"" % track_name)

EDIT_OPTS = (
    Opt("db_file", help="Track database file"),
    Opt("track_name", help="Current track name"),
    Opt(("-N", "--new-name"), default=None, help="New track name"),
    Opt(("-D", "--description"), default=None, help="New track description"),
    Opt(("-F", "--fast5-files"), "fast5_reader", type=comma_split),
    Opt(("-r", "--recursive"), "fast5_reader", action="store_true"),
)
def edit(conf, db=None):
    if db is None:
        db = TrackSQL(conf.db_file)
    track_id = db._verify_track(conf.track_name)

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
    db.cur.execute(query, params)

    if len(conf.fast5_files) > 0:
        _set_fast5s(track_id, conf.fast5_reader, db)
        
    db.con.commit()
    db.con.close()

def _set_fast5s(track_id, prms, db):
    fast5s = pd.read_sql_query(
        "SELECT fast5.id, filename FROM fast5 " \
        "JOIN read ON fast5.id = fast5_id " \
        "JOIN alignment ON read.id = read_id " \
        "WHERE track_id = ?",
        db.con, index_col="id", params=(track_id,))

    basenames = fast5s["filename"].map(os.path.basename)

    new_paths = {os.path.basename(path) : path for path in parse_fast5_paths(prms.fast5_files, prms.recursive)}

    fast5s["filename"] = basenames.map(new_paths)

    na = fast5s["filename"].isna()
    if na.any():
        missing = basenames[na].to_numpy() 
        if len(missing) > 10:
            missing = missing[:10] + ["..."]
        raise ValueError("Fast5 files found: " + ", ".join(missing))
    
    fast5s.to_sql(
        "fast5", con=db.con,
        if_exists="append",
        index_label="id")

MERGE_OPTS = (
    Opt("dbs", nargs="+", type=str, help="Database files to merge. Will write to the first file if \"-o\" is not specified. "),
    #Opt(("-n", "--track_names"), nargs="+", help="Names of tracks to merge. Will merge all tracks if not specified"),
    Opt(("-o", "--out-db"), type=str, default=None, help="Output database file. Will output to the first input file if not specified"),
)
def merge(conf):
    """Merge databases into a single file"""
    if conf.out_db is None:
        out_db = conf.dbs[0]
        in_dbs = conf.dbs[1:]
    else:
        in_dbs = conf.dbs
        out_db = conf.out_db

    db = TrackSQL(out_db)
    
    def max_id(table, field="id"):
        i = db.cur.execute(f"SELECT max({field}) FROM {table}").fetchone()[0]
        if i is None:
            return 0
        return i

    for fname in in_dbs:
        track_shift = max_id("track")
        aln_shift = max_id("alignment")+1
        fast5_shift = max_id("fast5")

        sys.stderr.write(f"Merging \"{fname}\"...\n")

        db.cur.execute(f"ATTACH DATABASE '{fname}' AS input")
        db.con.commit()
        db.cur = db.con.cursor()

        query = "INSERT INTO track "\
               f"SELECT id+{track_shift},name,desc,config "\
                "FROM input.track"
        #db._add_where(wheres,params,"name",conf.track_names)
        #query = db._join_query(query, wheres)
        db.cur.execute(query)

        query = "INSERT INTO fast5 "\
               f"SELECT id+{fast5_shift},filename FROM input.fast5"
        db.cur.execute(query)

        query = "INSERT OR IGNORE INTO read "\
               f"SELECT id,fast5_id+{fast5_shift} FROM input.read"
        db.cur.execute(query)

        query = "INSERT INTO alignment "\
               f"SELECT id+{aln_shift},track_id+{track_shift},read_id,ref_name,ref_start,ref_end,fwd,samp_start,samp_end,tags "\
                "FROM input.alignment"
        db.cur.execute(query)

        query = "INSERT INTO dtw "\
               f"SELECT mref,aln_id+{aln_shift},start,length,current,stdv "\
                "FROM input.dtw"
        db.cur.execute(query)

        query = "INSERT INTO bcaln "\
               f"SELECT mref,aln_id+{aln_shift},start,length,bp,error "\
                "FROM input.bcaln"
        db.cur.execute(query)

        query = "INSERT INTO cmp "\
               f"SELECT mref,aln_a+{aln_shift},aln_b+{aln_shift},group_b,mean_ref_dist,jaccard "\
                "FROM input.cmp"
        db.cur.execute(query)

        db.con.commit()
        db.cur.execute("DETACH input")
        db.con.commit()
    db.close()


        
SUBCMDS = [
    (ls, LS_OPTS), 
    (delete, DELETE_OPTS), 
    (merge, MERGE_OPTS), 
    (edit, EDIT_OPTS), 
]


