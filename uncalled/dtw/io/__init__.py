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


import _uncalled

from ..aln_track import AlnTrack
from ...fast5 import parse_fast5_paths
from ...pore_model import PoreModel
from ...signal_processor import ProcessedRead

INPUT_PARAMS = np.array(["db_in", "eventalign_in", "tombo_in"])
OUTPUT_PARAMS = np.array(["db_out", "eventalign_out"])

class TrackIO:
    def __init__(self, filename, conf, model, mode):
        self.conf = conf
        self.prms = self.conf.tracks.io
        self.model = model

        self.tracks = list()

        if not hasattr(self, "prev_aln_id"):
            self.prev_aln_id = 0

        if filename is None:
            self.filename = None
            self.track_names = None
        else:
            spl = filename.split(":")

            if len(spl) == 1:
                self.filename = filename
                self.track_names = None

            elif len(spl) == 2:
                self.filename = spl[0]
                self.track_names = spl[1].split(",")
            else:
                raise ValueError("Invalid database specifier format: " + filename)

        if mode == "w":
            self.write_mode = True
        elif mode == "r":
            self.write_mode = False
        else:
            raise ValueError("TrackIO mode must be either \'w\' or \'r\'")

    def init_alignment(self, read_id, fast5):
        if fast5 == self.prev_fast5[0]:
            fast5_id = self.prev_fast5[1]
        else:
            fast5_id = self.init_fast5(fast5)
            self.prev_fast5 = (fast5, fast5_id)

        if self.prev_read != read_id:
            self.init_read(read_id, fast5_id)
            self.prev_read = read_id

        self.prev_aln_id += 1
        return self.prev_aln_id

    def init_write_mode(self):
        if self.track_names is None:
            name = os.path.splitext(os.path.basename(self.filename))[0]
            self.track_names = [name]

        if len(self.track_names) > 1:
            raise ValueError("Can only write to one track")
        name = self.track_names[0]

        self.prev_fast5 = (None, None)
        self.prev_read = None

        track = AlnTrack(self, None, name, name, self.conf)
        self.tracks.append(track)

        return track

from .sqlite import TrackSQL
from .eventalign import Eventalign

class TrackHDF5(TrackIO):
    FORMAT = "hdf5"
    def __init__(self, filename, mode, conf):
        TrackIO.__init__(self, filename, mode, conf)

        new_file = not os.path.exists(filename)

        self.db = pd.HDFStore(filename)
        self.open = True


    def close(self):
        if self.open:
            self.db.close()
            self.open = False

    def init_tables(self):
        pass

    def init_write_mode(self):
        pass

    def init_track(self, track):
        pass

    def write_alignment(self, aln_df):
        pass 

    def init_fast5(self, filename):
        return fast5_id

    def init_read(self, read_id, fast5_id):
        pass 

    def write_layers(self, df, index=["pac","aln_id"], read=None):
        pass 

    def get_fast5_index(self, track_id=None):
        pass 

    def query_track(self, name=None):
        pass

    def query_alignments(self, track_id=None, read_id=None, aln_id=None, coords=None, full_overlap=False, order=None, chunksize=None):
        pass
        
    def query_compare(self, layers, track_id=None, coords=None, aln_id=None):
        pass

    def query_layers(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, order=["pac"], chunksize=None, full_overlap=False):
        pass

def _db_track_split(db_str):
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


_LS_QUERY = "SELECT name,desc,COUNT(alignment.id) FROM track " \
            "JOIN alignment ON track.id == track_id GROUP BY name"
def ls(conf, db=None):
    if db is None:
        db = TrackSQL(conf, None, "r")
    print("\t".join(["Name", "Description", "Alignments"]))

    for row in db.cur.execute(_LS_QUERY).fetchall():
        print("\t".join(map(str, row)))

    db.con.commit()

    #db.close()

def delete(track_name=None, db=None, conf=None):
    if db is None:
        db = TrackSQL(conf, None, "r")


    if track_name is None:
        track_name = conf.track_name

    db._verify_track(track_name)

    db.cur.execute("DELETE FROM track WHERE name == ?", (track_name,))
    db.con.commit()
    print("Deleted track \"%s\"" % track_name)

def edit(conf, db=None):
    fast5s = _uncalled._Fast5Reader.Params(conf.fast5_reader)
    fast5_change = len(conf.fast5_files) > 0
    track_name = conf.track_name
    if db is None:
        db = TrackSQL(conf, None, "r")
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
        conf.fast5_reader = fast5s
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

def _set_fast5s(track_id, fast5_files, db):
    fast5s = pd.read_sql_query(
        "SELECT fast5.id, filename FROM fast5 " \
        "JOIN read ON fast5.id = fast5_id " \
        "JOIN alignment ON read.id = read_id " \
        "WHERE track_id = ?",
        db.con, index_col="id", params=(track_id,))

    basenames = fast5s["filename"].map(os.path.basename)

    new_paths = {os.path.basename(path) : path for path in parse_fast5_paths(fast5_files, True)}

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

    db = TrackSQL(conf, None, "w")
    
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
               f"SELECT pac,aln_id+{aln_shift},start,length,current,stdv,kmer "\
                "FROM input.dtw"
        db.cur.execute(query)

        query = "INSERT INTO bcaln "\
               f"SELECT pac,aln_id+{aln_shift},start,length,bp,error "\
                "FROM input.bcaln"
        db.cur.execute(query)

        query = "INSERT INTO band "\
               f"SELECT aln_id+{aln_shift},pac,pac_end,sample_start,sample_end "\
                "FROM input.band"
        db.cur.execute(query)

        query = "INSERT INTO cmp "\
               f"SELECT pac,aln_a+{aln_shift},aln_b+{aln_shift},group_b,mean_ref_dist,jaccard "\
                "FROM input.cmp"
        db.cur.execute(query)

        db.con.commit()
        db.cur.execute("DETACH input")
        db.con.commit()
    db.close()




