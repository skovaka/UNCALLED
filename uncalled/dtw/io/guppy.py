import numpy as np
import pandas as pd
import sys, os, re
from collections import defaultdict
from ..aln_track import AlnTrack
from ...signal_processor import ProcessedRead
from ...index import str_to_coord, RefCoord
from ...pore_model import PoreModel
from ...fast5 import Fast5Reader
from ... import Config
from ont_fast5_api.fast5_interface import get_fast5_file
from . import TrackIO, OUT_EXT, OUTPUT_PARAMS
import _uncalled

bam_re = re.compile("bam_runid_([^_]+)_(\d+)_\d+\.bam")

class Guppy(TrackIO):
    FORMAT = "guppy"

    def __init__(self, filename, write, conf):
        self.dir = filename
        TrackIO.__init__(self, filename, write, conf)
        
        if not os.path.isdir(self.dir):
            raise FileNotFoundError(f"Guppy directory does not exist: {self.dir}")

        #if len(self.conf.fast5_reader.fast5_index) > 0:
        #    seqsum_fname = self.conf.fast5_reader.fast5_index
        #else:
        seqsum_fname = os.path.join(self.dir, "sequencing_summary.txt")
        if not os.path.isfile(seqsum_fname):
            raise FileNotFoundError(f"Guppy directory must contain sequencing_summary.txt")

        self.seqsum = pd.read_csv(seqsum_fname, sep="\t", usecols=["filename", "run_id", "batch_id"], chunksize=4000)

        self.has_fast5 = False
        self.bams = defaultdict(list)
        self.fast5_paths = dict()
        self.batches = dict()

        for root, dirs, files in os.walk(self.dir):
            for fname in files:
                if fname.endswith(".bam"):
                    m = bam_re.match(fname)
                    
                    #TODO need to support multiple BAM inputs to load failed reads
                    if "fail" in fname:
                         continue

                    if m is None:
                        raise RuntimeError(f"Unknown Guppy BAM filename format: {fname}")
                    runid = m.group(1)
                    batch = int(m.group(2))
                    self.bams[(runid,batch)].append(os.path.join(root, fname))
                elif fname.endswith(".fast5"):
                    self.fast5_paths[fname] = root

        out_prms = [getattr(self.conf.tracks.io, p) is not None for p in OUTPUT_PARAMS]
        if np.sum(out_prms) != 1:
            raise ValueError("Must specify one output format")
        self.out_fmt = OUTPUT_PARAMS[out_prms][0]
        self.out_prefix = getattr(self.conf.tracks.io, self.out_fmt)
        if not os.path.isdir(self.out_prefix):
            self.out_prefix += "."

        #fast5_dir = os.path.join(self.dir, "workspace")
        #if self.has_fast5 and os.path.isdir(fast5_dir):
        #    self.conf.fast5_reader.fast5_files = [fast5_dir]
        #    self.conf.fast5_reader.recursive = True

        self.fast5_in = None
        self.read_id_in = None

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

    def init_write_mode(self):
        raise RuntimeError("Writing in guppy format not supported")

    def init_read_mode(self):
        name = self.filename

        self.conf.fast5_reader.fast5_files = [self.prms.guppy_in]

        self.init_track(1, name, name, self.conf.to_toml())

    def iter_batches(self):
        self.batches = set()
        for chunk in self.seqsum:
            for _,row in chunk.iterrows():
                batch = (row["run_id"], row["batch_id"])
                fast5_dir = self.fast5_paths[row["filename"]]
                if not batch in self.batches:
                    self.batches.add(batch)

                    conf = Config(self.conf)
                    conf.tracks.io.guppy_in = None
                    conf.tracks.io.bam_in = self.bams[batch][0]
                    conf.fast5_files = [os.path.join(fast5_dir, row["filename"])]

                    sbatch = ".".join(map(str, batch))
                    setattr(conf.tracks.io, self.out_fmt, f"{self.out_prefix}{sbatch}.{OUT_EXT[self.out_fmt]}")

                    yield conf


    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
            yield alns,layers
            
