"""Edit, merge, and ls alignment databases
        print(d)

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
from collections import namedtuple
import multiprocessing as mp
import pysam
from time import time

from ..aln_track import AlnTrack
from ...config import Config
from ... import RefCoord, ExceptionWrapper
from ...read_index import ReadIndex

INPUT_PARAMS = np.array(["sql_in", "eventalign_in", "tombo_in", "bam_in"])
OUTPUT_PARAMS = np.array(["sql_out", "tsv_out", "eventalign_out", "bam_out", "model_dir"])

OUT_EXT = {
    "sql_out" : "db", 
    "tsv_out" : "tsv", 
    "eventalign_out" : "txt", 
    "bam_out" : "bam"
}

#AlnTrack = namedtuple("AlnTrack", ["id", "name", "desc", "conf"])

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

        if not hasattr(self, "aln.id"):
            self.aln_id = 0

        if filename is None:
            self.filename = None
        elif isinstance(filename, str):
            self.filename = filename

        if self.prms.in_names is not None:
            self.in_tracks = self.prms.in_names
        else:
            self.in_tracks = None

        self.write_mode = write

    @property
    def read_filter(self):
        return self.tracks.read_index.read_filter

    def init_write_mode(self):
        if self.prms.out_name is not None:
            track_name = self.prms.out_name
        else:
            track_name = os.path.splitext(os.path.basename(self.filename))[0]

        self.prev_fast5 = (None, None)
        self.prev_read = None

        self.track_out = self.init_track(track_name, track_name, self.conf)

        self.out_buffer = None

    def init_track(self, name, desc, conf, id=None):
        if id is None:
            id = self.next_id
        self.next_id = id + 1

        t = AlnTrack(id, name, desc, conf)
        self.aln_tracks.append(t)
        self.conf.load_config(conf)

        return t

    def fill_tracks(self, coords, alignments, layers):
        layers = layers.droplevel(0)

        for track in self.aln_tracks:
            track_alns = alignments[alignments["track.name"] == track.name]
            i = layers.index.get_level_values("aln.id").isin(track_alns.index)
            track_layers = layers.iloc[i]

            track.set_data(coords, track_alns, track_layers)

    def _init_output(self, buffered, mode="w"):
        if buffered:
            self.output = None
            self.out_buffer = list()
        else:
            if self.filename == "-":
                self.output = sys.stdout
            else:
                self.output = open(self.filename, mode)

    def _set_output(self, out):
        if self.prms.buffered:
            self.out_buffer.append(out)
        else:
            self.output.write(out)

    def write_buffer(self, buf=None):
        if buf is None:
            buf = [self.out_buffer]

        for out in buf:
            self.output.write(out)
            
    def next_aln_id(self):
        ret = self.aln_id
        self.aln_id += 1
        return ret

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
from .model_trainer import ModelTrainer

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
    "model_dir" : ModelTrainer,
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
    conf.tracks.layers = ["dtw"]
    
    if conf.tracks.io.tombo_in is not None:
        conf.read_index.paths = conf.tracks.io.tombo_in

    if conf.tracks.io.processes == 1:
        convert_single(conf)
    elif conf.tracks.io.bam_in is not None:
        convert_pool(conf)
    else:
        raise ValueError("Parallel convert is only supported for BAM input") 


def _init_tracks(conf):
    from .. import Tracks

    tracks = Tracks(conf=conf)

    if len(tracks.inputs) == 2:
        if tracks.bam_in is None:
            raise ValueError("Only one non-BAM input can be specified")
        ignore_bam = True
    else:
        ignore_bam = False

    return tracks, ignore_bam

def convert_pool(conf):
    tracks,ignore_bam = _init_tracks(conf)

    def iter_args(): 
        i = 0
        aln_count = 0
        for read_ids, sam_strs in tracks.bam_in.iter_str_chunks():
            reads = ReadIndex()#tracks.read_index.subset(read_ids)
            yield (tracks.conf, tracks.model, sam_strs, reads, aln_count, tracks.bam_in.header)
            aln_count += len(sam_strs)

    try:
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=conf.tracks.io.processes, maxtasksperchild=4) as pool:

            i = 0
            if conf.tracks.io.ordered_out:
                itr = pool.imap
            else:
                itr = pool.imap_unordered

            for out in itr(convert_worker, iter_args(), chunksize=1):
                tracks.output.write_buffer(out)

    except Exception as e:
        raise ExceptionWrapper(e).re_raise()

def convert_worker(args):
    from .. import Tracks
    conf,model,sam_strs,reads,aln_start,header = args

    conf.tracks.io.buffered = True
    #conf.tracks.io.bam_in = None

    header = pysam.AlignmentHeader.from_dict(header)

    tracks = Tracks(model=model, read_index=reads, conf=conf)
    for s in sam_strs:
        sam = pysam.AlignedSegment.fromstring(s, header)
        aln = tracks.bam_in.sam_to_aln(sam, load_moves=False)
        aln.instance.id += aln_start
        tracks.write_alignment(aln)

    tracks.close()
    return tracks.output.out_buffer
        

def convert_single(conf):
    tracks, ignore_bam = _init_tracks(conf)
    for aln in tracks.iter_reads(ignore_bam=ignore_bam):
        tracks.write_alignment(aln)

    tracks.close()

