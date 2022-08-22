import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack, LAYER_META, parse_layers
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

        print(self.prms.tsv_cols)
        self.columns = pd.MultiIndex.from_tuples(parse_layers(self.prms.tsv_cols))
        self._header = True

    def init_read_mode(self):
        raise RuntimeError("Reading from TSV not yet supported")
        #name = self.track_names[0]
        #t = AlnTrack(self, None, name, name, self.conf, self.model)
        #self.tracks.append(t)

    def write_layers(self, track):
        df = track.layers_desc_index#.reset_index()
        df = df[self.columns.intersection(df.columns)].dropna(how="all", axis=0).reset_index()

        if self._header:
            self.columns = df.columns
            labels = [
                ".".join([c for c in col if len(c) > 0]) 
                for col in self.columns]
            self.out.write("\t".join(labels) + "\n")
            self._header = False

        df = df.loc[:,self.columns]

        self.out.write(df.to_csv(sep="\t", header=False, na_rep=self.prms.tsv_na, index=False))

    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        self.out.close()

