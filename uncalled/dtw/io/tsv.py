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

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

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

        TrackIO._init_output(self, self.prms.buffered)

    def init_read_mode(self):
        raise RuntimeError("Reading from TSV not yet supported")

    def write_layers(self, track, groups):
        track.calc_layers(self.columns)
        df = track.layers_desc_index#.reset_index()
        
        df = df[self.columns.intersection(df.columns)].dropna(how="all", axis=0)

        for col in df.columns[df.columns.get_level_values(-1).str.endswith("kmer")]:
            kmers = df[col].dropna()
            df[col] = track.model.kmer_to_str(kmers)


        for col in self.extras:
            df[col] = EXTRA_FNS[col](track,df)

        df.reset_index(inplace=True, drop=self.prms.tsv_noref)

        if self._header:
            self.columns = df.columns
            labels = [
                ".".join([c for c in col if len(c) > 0]) 
                for col in self.columns]
            self._header = False
            if not self.prms.buffered:
                self.output.write("\t".join(labels) + "\n")

        df = df.loc[:,self.columns]

        out = df.to_csv(sep="\t", header=False, na_rep=self.prms.tsv_na, index=False)

        self._set_output(out)


    def close(self):
        if not self.prms.buffered:
            self.output.close()

