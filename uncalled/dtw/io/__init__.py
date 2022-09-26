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

INPUT_PARAMS = np.array(["sql_in", "eventalign_in", "tombo_in", "bam_in"])
OUTPUT_PARAMS = np.array(["sql_out", "tsv_out", "eventalign_out", "bam_out"])

OUT_EXT = {
    "sql_out" : "db", 
    "tsv_out" : "tsv", 
    "eventalign_out" : "txt", 
    "bam_out" : "bam"
}

class TrackIO:
    def __init__(self, filename, conf, mode):
        self.conf = conf
        self.prms = self.conf.tracks.io

        self.read = None

        self.tracks = None #list()

        if not hasattr(self, "prev_aln_id"):
            self.prev_aln_id = 0

        if filename is None:
            self.filename = None
            self.track_names = None
        elif isinstance(filename, str):
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

    def init_alignment(self, read_id, fast5, read):
        self.read = read
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

        return self.init_track(None, name, name, self.conf.to_toml())

    def init_track(self, id, name, desc, conf):
        row = pd.DataFrame({
            "id" : id,
            "name" : name,
            "desc" : desc,
            "config" : conf
        }, index=[0])
        self.init_tracks(row)
        return row

    def init_tracks(self, df):
        if self.tracks is None:
            self.tracks = df
        else:
            self.tracks = pd.concat([self.tracks, df])
            
    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def query_compare(self, layers, track_id=None, coords=None, aln_id=None):
        pass



from .sqlite import TrackSQL
from .tsv import TSV
from .bam import BAM
from .eventalign import Eventalign
from .tombo import Tombo
from .guppy import Guppy

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


def convert(conf):
    """Convert between signal alignment file formats"""
    from .. import Tracks

    conf.tracks.layers = ["dtw"]
    tracks = Tracks(conf=conf)

    for read_id, read in tracks.iter_reads():
        sys.stderr.write(f"{read_id}\n")
        aln = read.alns[0]
        read.init_alignment(read_id, read.get_read_fast5(read_id), aln.coords, {"dtw" : aln.layers["dtw"].droplevel(1)})
        
        read.write_alignment()

    tracks.close()

