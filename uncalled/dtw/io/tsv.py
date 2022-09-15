import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from . import TrackIO
import _uncalled

EXTRA_FNS = {
    "read_id" : lambda track,df: track.alignments.loc[df.index.get_level_values("aln_id"), "read_id"].to_numpy()
}

EXTRA_COLS = pd.Index(EXTRA_FNS.keys())

class TSV(TrackIO):
    FORMAT = "tsv"

    def __init__(self, conf, mode):
        filename = conf.tracks.io.tsv_out
        TrackIO.__init__(self, filename, conf, mode)

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

        layer_cols = pd.Index(self.prms.tsv_cols).difference(EXTRA_COLS)
        self.extras = EXTRA_COLS.intersection(self.prms.tsv_cols)

        layers = list(parse_layers(layer_cols, add_deps=False))
        self.columns = pd.MultiIndex.from_tuples(layers)
        self.conf.tracks.layers += layers
        self._header = True

    def init_read_mode(self):
        raise RuntimeError("Reading from TSV not yet supported")
        #name = self.track_names[0]
        #t = AlnTrack(self, None, name, name, self.conf, self.model)
        #self.tracks.append(t)

    def write_layers(self, track, groups):
        track.calc_layers(self.columns)
        df = track.layers_desc_index#.reset_index()
        
        df = df[self.columns.intersection(df.columns)].dropna(how="all", axis=0)

        for col in df.columns[df.columns.get_level_values(-1).str.endswith("kmer")]:
            kmers = df[col].dropna()
            df.loc[kmers.index, col] = track.model.kmer_to_str(kmers)


        for col in self.extras:
            df[col] = EXTRA_FNS[col](track,df)
        df.reset_index(inplace=True)

        if self._header:
            self.columns = df.columns
            labels = [
                ".".join([c for c in col if len(c) > 0]) 
                for col in self.columns]
            self.out.write("\t".join(labels) + "\n")
            self._header = False

        df = df.loc[:,self.columns]

        self.out.write(df.to_csv(sep="\t", header=False, na_rep=self.prms.tsv_na, index=False))

    def close(self):
        self.out.close()

