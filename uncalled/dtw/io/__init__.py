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
OUTPUT_PARAMS = np.array(["db_out", "tsv_out", "eventalign_out"])

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
        #track = AlnTrack(self, None, name, name, self.conf)
        #self.tracks.append(track)
        #return track

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
            
        


from .sqlite import TrackSQL
from .tsv import TSV
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



