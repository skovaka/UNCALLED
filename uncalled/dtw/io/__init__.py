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

from ..aln_track import AlnTrack
from ...config import Config

INPUT_PARAMS = np.array(["sql_in", "eventalign_in", "tombo_in", "bam_in"])
OUTPUT_PARAMS = np.array(["sql_out", "tsv_out", "eventalign_out", "bam_out"])

OUT_EXT = {
    "sql_out" : "db", 
    "tsv_out" : "tsv", 
    "eventalign_out" : "txt", 
    "bam_out" : "bam"
}

class TrackIO:
    def __init__(self, filename, write, tracks, track_count):
        if isinstance(tracks, Config):
            self.tracks = None
            self.conf = tracks
        else:
            self.tracks = tracks
            self.conf = tracks.conf
        self.prms = self.conf.tracks.io
        self.next_id = track_count+1

        self.read = None
        self.bam = None

        self.aln_tracks = list()

        if not hasattr(self, "prev_aln_id"):
            self.prev_aln_id = 1

        if filename is None:
            self.filename = None
        elif isinstance(filename, str):
            self.filename = filename

        if self.prms.in_names is not None:
            self.in_tracks = self.prms.in_names
        else:
            self.in_tracks = None

            #spl = filename.split(":")

            #if len(spl) == 1:
            #    self.filename = filename
            #    self.track_names = None

            #elif len(spl) == 2:
            #    self.filename = spl[0]
            #    self.track_names = spl[1].split(",")
            #else:
            #    raise ValueError("Invalid database specifier format: " + filename)

        #self.track_names = None
        self.write_mode = write

    def init_alignment(self, read_id, fast5, read, bam):
        self.read = read
        self.bam = bam

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
        if self.prms.out_name is not None:
            self.out_track = self.prms.out_name
        else:
            self.out_track = os.path.splitext(os.path.basename(self.filename))[0]

        self.prev_fast5 = (None, None)
        self.prev_read = None

        self.out_id = self.init_track(self.out_track, self.out_track, self.conf)

    def init_track(self, name, desc, conf, id=None):
        if id is None:
            id = self.next_id
        self.next_id = id + 1

        self.aln_tracks.append(AlnTrack(id, name, desc, conf))
        self.conf.load_config(conf)

        return id

    def fill_tracks(self, coords, alignments, layers):
        layers = layers.droplevel(0)

        for track in self.aln_tracks:
            track_alns = alignments[alignments["track_id"] == track.id]
            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i]

            track.set_data(coords, track_alns, track_layers)
            
    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def query_compare(self, layers, track_id=None, coords=None, aln_id=None):
        pass



from .sqlite import TrackSQL as SQL
from .tsv import TSV
from .bam import BAM
from .eventalign import Eventalign
from .tombo import Tombo
from .guppy import Guppy

INPUTS = {
    "sql_in" : SQL, 
    "bam_in" : BAM,
    "eventalign_in" : Eventalign, 
    "tombo_in" : Tombo, 
}

OUTPUTS = {
    "sql_out" : SQL,
    "bam_out" : BAM,
    "eventalign_out" : Eventalign,
    "tsv_out" : TSV,
}

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
    conf.tracks.load_fast5s = True
    tracks = Tracks(conf=conf)

    for read_id, read in tracks.iter_reads():
        sys.stderr.write(f"{read_id}\n")
        aln = read.alns[0]
        dtw = aln.layers["dtw"].droplevel(1)
        read.init_alignment(read_id, read.get_read_fast5(read_id), aln.coords, {"dtw" : dtw})
        
        read.write_alignment()

    tracks.close()

