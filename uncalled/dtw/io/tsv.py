import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from . import TrackIO
import _uncalled

class TSV(TrackIO):
    FORMAT = "tsv"

    def __init__(self, conf, model, mode):
        filename = conf.tracks.io.tsv_out
        TrackIO.__init__(self, filename, conf, model, mode)

        if self.filename == "-":
            self.out = sys.stdout
        else:
            self.out = open(self.filename, mode)

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

    def init_write_mode(self):
        TrackIO.init_write_mode(self)
        
        self._header = True

    def init_read_mode(self):
        raise RuntimeError("Reading from TSV not yet supported")
        #name = self.track_names[0]
        #t = AlnTrack(self, None, name, name, self.conf, self.model)
        #self.tracks.append(t)

    def write_layers(self, df):
        df = df.reset_index()
        if self._header:
            columns = [
                ".".join([c for c in col if len(c) > 0]) 
                for col in df.columns]
            self.out.write("\t".join(columns) + "\n")
            self._header = False
        self.out.write(df.to_csv(sep="\t", header=False, na_rep=self.prms.tsv_na, index=False))

    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        self.out.close()

