import sqlite3
import os
import collections
import numpy as np
import pandas as pd
import sys

from ..aln_track import AlnTrack
from ..layers import LAYER_META, LAYER_DB_GROUPS
from ... import config
from . import TrackIO

class TrackSQL(TrackIO):
    FORMAT = "db"
    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

        self.track_shift = track_count#+1

        new_file = not os.path.exists(self.filename)
        self.con = sqlite3.connect(self.filename)
        self.cur = self.con.cursor()
        self.cur.execute("PRAGMA foreign_keys = ON")
        self.con.commit()
        self.open = True

        #if new_file:
        self.init_tables()

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

        prev_id = self.cur.execute("SELECT MAX(id) FROM alignment").fetchone()[0]
        if prev_id is not None:
            self.prev_aln_id = prev_id

    def close(self):
        if self.open:
            self.cur.execute("CREATE INDEX IF NOT EXISTS aln_read_idx ON alignment (read_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS aln_fwd_idx ON alignment (fwd);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS dtw_idx ON dtw (pac, aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS dtw_aln_idx ON dtw (aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS moves_idx ON moves (pac, aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS moves_aln_idx ON moves (aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS band_idx ON band (pac, aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS band_aln_idx ON band (aln_id);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS cmp_idx ON cmp (pac, aln_a, aln_b, group_b);")
            self.cur.execute("CREATE INDEX IF NOT EXISTS cmp_aln_idx ON cmp (aln_a);")
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
                pac INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                current REAL,
                stdv REAL,
                events REAL DEFAULT 1,
                kmer INTEGER,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS moves (
                pac INTEGER,
                aln_id INTEGER,
                start INTEGER,
                length INTEGER,
                indel INTEGER,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
                #kmer INTEGER,
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS cmp (
                pac INTEGER,
                aln_a INTEGER,
                aln_b INTEGER,
                group_b TEXT DEFAULT "dtw",
                dist REAL, 
                jaccard REAL, 
                FOREIGN KEY (aln_a) REFERENCES alignment (id) ON DELETE CASCADE,
                FOREIGN KEY (aln_b) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS band (
                aln_id INTEGER,
                pac INTEGER,
                pac_end INTEGER,
                sample_start INTEGER,
                sample_end INTEGER,
                FOREIGN KEY (aln_id) REFERENCES alignment (id) ON DELETE CASCADE
            );""")
        #self.con.commit()

    def init_write_mode(self):
        TrackIO.init_write_mode(self)

        for table in ["dtw", "moves", "cmp", "band"]:
            self.cur.execute("DROP INDEX IF EXISTS %s_idx" % table)

        track = self.aln_tracks[0]

        if self.prms.init_track:
            try:
                tid = self.write_track(track)
            except Exception as err:
                if len(self.query_track(track.name)) > 0:
                    if self.prms.append:
                        pass
                    elif self.prms.overwrite:
                        sys.stderr.write("Deleting existing track...\n")
                        delete(track.name, self)
                        tid = self.write_track(track)
                    else:
                        raise ValueError(f"Database already contains track named \"{ track.name}\". Specify a different name, write to a different file")
                else:
                    raise err
        else:
            tid = 0

        self.write_id = tid

        #self.aln_tracks[0].id = tid

    def init_read_mode(self):
        self.table_columns = dict()
        for table in LAYER_DB_GROUPS:
            info = pd.read_sql_query(f"PRAGMA table_info({table})", self.con)
            if not table in LAYER_META.index.get_level_values(0): continue
            layers = LAYER_META.loc[table].query("base").index
            missing = layers.difference(info["name"])
            if len(missing) > 0:
                n = len(missing)
                missing = "\", \"".join(missing)
                sys.stderr.write(f"Warning: \"{table}\" table missing {n} columns (\"{missing}\")\n")
            self.table_columns[table] = layers.intersection(info["name"])

        df = self.query_track(self.prms.in_names).set_index("name").reindex(self.prms.in_names)

        missing = df["id"].isnull()
        if missing.any():
            bad_names = df[missing].index
            all_names = self.query_track()["name"]
            raise ValueError("Alignment track not found: \"%s\" (tracks in database: \"%s\")" %
                             ("\", \"".join(bad_names), "\", \"".join(all_names)))

        #self.init_tracks(df.reset_index())

        for i,row in df.reset_index().iterrows():
            if "config" in row.index:
                conf = config.Config(toml=row["config"])
            else:
                conf = self.conf
            self.init_track(row["name"], row["desc"], conf, row["id"]+self.track_shift)

        #if self.tracks is None:
        #    self.tracks = df
        #else:
        #    self.tracks = pd.concat([self.tracks, df])


    def write_track(self, track):
        #self.cur.execute(
        #    "INSERT INTO track (id,name,desc,config) VALUES (?,?,?,?)",
        #    (track.id, track.name, track.desc, track.conf.to_toml())
        #)
        self.cur.execute(
            "INSERT INTO track (name,desc,config) VALUES (?,?,?)",
            (track.name, track.desc, track.conf.to_toml())
        )
        track_id = self.cur.lastrowid
        return track_id

    def update_config(self, track):
        self.cur.execute(
            "UPDATE track SET config = ? WHERE name = ?",
            (track.conf.to_toml(), track.name)
        )


    def write_alignment(self, aln_df):
        aln_df["track.id"] = self.write_id
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

    def write_layers(self, track, groups):
        layers = track.layers_pac_index
        for group in groups:
            df = layers[group].dropna(axis=0, how="all")
            base_layers = df.columns.intersection(LAYER_META[LAYER_META["base"]].loc[group].index)
            self.write_layer_group(group, df[base_layers])

    def write_layer_group(self, group, df):
        if group == "mvcmp": group = "cmp"
        if group == "cmp":
            df = df.droplevel("aln.id")
        df.to_sql(
            group, self.con, 
            if_exists="append", 
            method="multi", chunksize=999//len(df.columns),
            index=True)

    def get_fast5_index(self, track_id=None):
        query = "SELECT read.id AS read_id, filename FROM read " +\
                "JOIN fast5 ON fast5.id = fast5_id"

        wheres = list()
        params = list()

        if track_id is not None:
            query += " JOIN alignment ON read.id = alignment.read_id"
            self._add_where(wheres, params, "track.id", track_id)
        
        query = self._join_query(query, wheres)
        return pd.read_sql_query(query, self.con, params=params)

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
        select = "SELECT alignment.id, track_id, read_id, ref_name, ref_start, ref_end, fwd, samp_start, samp_end, tags, filename AS fast5 FROM alignment JOIN read ON read_id == read.id JOIN fast5 ON fast5.id = fast5_id"
        wheres = list()
        params = list()

        self._add_where(wheres, params, "alignment.id", aln_id)
        self._add_where(wheres, params, "track.id", track_id)
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

        df = pd.read_sql_query(query, self.con, params=params, chunksize=chunksize)
        df["track.id"] += self.track_shift #self.next_id+1
        return df.set_index(["track.id","id"])
        
    def _join_query(self, select, wheres, order=None):
        query = select
        if len(wheres) > 0:
            query += " WHERE " + " AND ".join(wheres)
        if order is not None:
            query += " ORDER BY " + ", ".join(order)
        return query

    def query_compare(self, layers, track_id=None, coords=None, aln_id=None):
        dtw = False
        moves = False
        fields = {"seq.pac", "aln_a", "aln_b", "group_b"}
        for group,layer in layers:
            if group == "cmp": dtw = True
            if group == "mvcmp": moves = True
            fields.add(layer)

        fields = ", ".join(fields)
        select = f"SELECT {fields} FROM cmp"

        wheres = list()
        params = list()

        if track_id is not None:
            select += " JOIN alignment ON id = cmp.aln_a"
            self._add_where(wheres, params, "track.id", track_id)

        if dtw and not moves:
            self._add_where(wheres, params, "group_b", "dtw")
        elif moves and not dtw:
            self._add_where(wheres, params, "group_b", "moves")

        self._add_where(wheres, params, "aln_a", aln_id)

        if coords is not None:
            wheres.append("(pac >= ? AND pac <= ?)")
            params += [str(coords.pacs.min()), str(coords.pacs.max())]

        query = self._join_query(select, wheres)
        return pd.read_sql_query(
            query, self.con, 
            index_col=["seq.pac", "aln_a", "aln_b", "group_b"], 
            params=params)

    def iter_refs(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, chunksize=None, full_overlap=False):
        if fwd is not None:
            strands = [int(fwd)]
        else:
            strands = [1, 0]

        for fwd in strands:
            itr = self.query_layers(layers, track_id, coords, aln_id, read_id, fwd, ["seq.pac"], chunksize, full_overlap)

            ret_alns = pd.DataFrame()
            chunks = list()
            prev_pac = None

            for layers in itr:
                layers.sort_index(inplace=True)
                next_pac = layers.index.get_level_values("seq.pac")[-1]
                #if next_pac 

                if prev_pac is None or prev_pac == next_pac:
                    chunks.append(layers)
                else:
                    chunks.append(layers.loc[slice(None),fwd,:next_pac-1])
                    ret_layers = pd.concat(chunks)
                    chunks = [layers.loc[slice(None),fwd,next_pac:]]

                    aln_ids = ret_layers.index.get_level_values("aln.id").unique()
                    if len(ret_alns) > 0:
                        aln_ids = aln_ids.difference(ret_alns.index)

                    ret_alns = self.query_alignments(track_id, aln_id=list(aln_ids))

                    yield ret_alns, ret_layers

                prev_pac = next_pac

            if len(chunks) > 0:
                ret_layers = pd.concat(chunks)
                chunks = [layers.loc[slice(None),fwd,next_pac]]

                aln_ids = ret_layers.index.get_level_values("aln.id").unique()
                if len(ret_alns) > 0:
                    aln_ids = aln_ids.difference(ret_alns.index)

                ret_alns = self.query_alignments(track_id, aln_id=list(aln_ids))

                yield ret_alns, ret_layers

    def query(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, order=["read_id", "seq.pac"], full_overlap=False):
        
        if track_id is not None:
            track_id = np.array(track_id)-self.track_shift

        layers = self.query_layers(layers, track_id, coords, aln_id, read_id, fwd, order, None, full_overlap)

        aln_ids = layers.index.unique("aln.id").to_numpy()
        alignments = self.query_alignments(aln_id=aln_ids)

        #self.fill_tracks(coords, alignments, layers)
        return alignments, layers

    def query_layers(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, order=["read_id", "seq.pac"], chunksize=None, full_overlap=False):

        group_layers = collections.defaultdict(list)
        renames = dict()
        fields = list()
        tables = list()
        for group,layer in layers:
            if not layer in self.table_columns[group]: 
                continue
            if group == "mvcmp":
                group = "cmp"
            field = group + "." + layer
            name = group + "_" + layer
            group_layers[group].append(name)
            renames[name] = layer
            fields.append(field + " AS " + name)
            if group not in tables:
                tables.append(group)

        fields += ["%s.%s AS idx_%s" % (tables[0],idx,idx) for idx in ("seq.pac", "aln.id")]
        fields.append("seq.fwd")
        fields.append("track.id")

        select = "SELECT " + ", ".join(fields) + " FROM " + tables[0]
        for table in tables[1:]:
            select += " LEFT JOIN %s ON %s.aln_id == idx_aln_id AND %s.pac == idx_pac" % ((table,)*3)

        wheres = list()
        params = list()

        #if track_id is not None or read_id is not None or "seq.fwd" in order:
        select += " JOIN alignment ON id = idx_aln_id"
        if track_id is not None:
            self._add_where(wheres, params, "track.id", track_id)
        if read_id is not None:
            self._add_where(wheres, params, "read_id", read_id)
        if fwd is not None:
            self._add_where(wheres, params, "seq.fwd", fwd)

        if "cmp" in group_layers:
            self._add_where(wheres, params, "group_b", "moves")

        if coords is not None:
            wheres.append("(idx_pac >= ? AND idx_pac <= ?)")
            params += [str(coords.pacs.min()), str(coords.pacs.max())]

            if full_overlap:
                wheres.append("ref_name = ? AND ref_start < ? AND ref_end > ?")
                params.append(coords.ref_name)
                params += [coords.refs.min(), coords.refs.max()]

        self._add_where(wheres, params, "idx_aln_id", aln_id)

        query = self._join_query(select, wheres, ["idx_"+o if o in {"aln.id","seq.pac"} else o for o in order])

        ret = pd.read_sql_query(
            query, self.con, 
            #index_col=["track.id", "seq.fwd", "idx_pac", "idx_aln_id"], 
            params=params, chunksize=chunksize)

        def make_groups(df):
            grouped = dict()
            df["track.id"] += self.track_shift #self.next_id+1
            df.set_index(["track.id", "seq.fwd", "idx_pac", "idx_aln_id"], inplace=True)
            df.index.names = ("track.id", "seq.fwd", "seq.pac", "aln.id")
            for group, layers in group_layers.items():
                gdf = df[layers]
                grouped[group] = df[layers].rename(columns=renames)
            df = pd.concat(grouped, names=("group", "layer"), axis=1)
            return df#.sort_index()
                
        if chunksize is None:
            return make_groups(ret)
        
        return (make_groups(df) for df in ret if len(df) > 0)

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
        layer_iter = self.query_layers(
            layers, track_id, coords, aln_id, read_id, fwd,
            ["read_id", "seq.pac"], self.prms.ref_chunksize, full_overlap) 

        layer_leftovers = pd.DataFrame()

        for layers in layer_iter:
            layers = pd.concat([layer_leftovers, layers])

            ids = layers.index.unique("aln.id").to_numpy()
            alignments = self.query_alignments(aln_id=ids, order=["read_id"])

            #print("HERE")
            last_read = alignments.iloc[-1]["read_id"]
            mask = alignments["read_id"] == last_read
            last_alns = alignments.index[mask].get_level_values("id")
            alignments = alignments[~mask]

            #print(last_alns)
            #print(layers)

            layer_leftovers = layers.loc[layers.index.get_level_values("aln.id").isin(last_alns)]
            #layer_leftovers = layers.loc[last_layers]

            layers = layers.drop(index=layer_leftovers.index)

            yield alignments, layers

        if len(layer_leftovers) > 0:
            ids = layer_leftovers.index.unique("aln.id").to_numpy()
            alignments = self.query_alignments(aln_id=ids, order=["read_id"])
            yield alignments, layer_leftovers
    
    def _verify_track(self, track_name):
        ids = self.cur.execute("SELECT id FROM track WHERE name == ?", (track_name,)).fetchall()
        if len(ids) == 0:
            raise ValueError(f"Track does not exist: \"{track_name}\"\n")
        return ids[0][0]

    def delete_track(self, track_name):
        track = self._verify_track(track_name)
        self.cur.execute("DELETE FROM track WHERE id == ?", (track,))
        self.con.commit()

_LS_QUERY = "SELECT name,desc,COUNT(alignment.id) FROM track " \
            "JOIN alignment ON track.id == track_id GROUP BY name"
def ls(conf, db=None):
    if db is None:
        db = TrackSQL(conf.tracks.io.sql_in, False, conf, 1)
    print("\t".join(["Name", "Description", "Alignments"]))

    for row in db.cur.execute(_LS_QUERY).fetchall():
        print("\t".join(map(str, row)))

    db.con.commit()

def delete(track_name=None, db=None, conf=None):
    if db is None:
        #db = TrackSQL(conf, "r")
        db = TrackSQL(conf.tracks.io.sql_in, False, conf)

    if track_name is None:
        track_name = conf.track_name

    db.delete_track(track_name)

    print("Deleted track \"%s\"" % track_name)

def edit(conf, db=None):
    fast5_change = len(conf.paths) > 0
    track_name = conf.track_name
    if db is None:
        #db = TrackSQL(conf, "r")
        db = TrackSQL(conf.tracks.io.sql_in, False, conf, 1)
    track_id = db._verify_track(track_name)

    updates = []
    params = []
    if conf.new_name:
        updates.append("name = ?")
        params.append(conf.new_name)
    if conf.description:
        updates.append("desc = ?")
        params.append(conf.description)
    if fast5_change:
        updates.append("config = ?")
        conf.read_index = fast5s
        del conf.track_name
        del conf.new_name
        del conf.description
        params.append(conf.to_toml())
    if len(updates) == 0:
        sys.stderr.write("No changes made. Must specify new name (-N) or description (-D)\n")
        return

    params.append(track_name)
        
    query = "UPDATE track SET " + ", ".join(updates) + " WHERE name == ?"
    db.cur.execute(query, params)

    if fast5_change:
        _set_fast5s(track_id, fast5s, db)
        
    db.con.commit()
    db.con.close()

def _set_fast5s(track_id, paths, db):
    fast5s = pd.read_sql_query(
        "SELECT fast5.id, filename FROM fast5 " \
        "JOIN read ON fast5.id = fast5_id " \
        "JOIN alignment ON read.id = read_id " \
        "WHERE track_id = ?",
        db.con, index_col="id", params=(track_id,))

    basenames = fast5s["filename"].map(os.path.basename)

    new_paths = {os.path.basename(path) : path for path in parse_fast5_paths(paths, True)}

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

def merge(conf):
    """Merge databases into a single file"""
    #if conf.out_db is None:
    #    out_db = conf.dbs[0]
    #    in_dbs = conf.dbs[1:]
    #else:
    in_dbs = conf.dbs
    #    out_db = conf.out_db

    conf.tracks.io.init_track = False

    #db = TrackSQL(conf, "w")
    db = TrackSQL(conf.tracks.io.sql_out, True, conf)
    
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
               f"SELECT pac,aln_id+{aln_shift},start,length,current,stdv,events,kmer "\
                "FROM input.dtw"
        db.cur.execute(query)

        query = "INSERT INTO moves "\
               f"SELECT pac,aln_id+{aln_shift},start,length,indel "\
                "FROM input.moves"
        db.cur.execute(query)

        query = "INSERT INTO band "\
               f"SELECT aln_id+{aln_shift},pac,pac_end,sample_start,sample_end "\
                "FROM input.band"
        db.cur.execute(query)

        query = "INSERT INTO cmp "\
               f"SELECT pac,aln_a+{aln_shift},aln_b+{aln_shift},group_b,dist,jaccard "\
                "FROM input.cmp"
        db.cur.execute(query)

        db.con.commit()
        db.cur.execute("DETACH input")
        db.con.commit()
    db.close()


