import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ...aln import LAYER_META, parse_layers
from . import TrackIO
import _uncalled

EXTRA_FNS = {
    "read_id" : lambda aln: aln.read_id
}

EXTRA_COLS = pd.Index(EXTRA_FNS.keys())
INDEX_COLS = [("seq","name"),("seq","pos"),("seq","strand"),("aln","id")]

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

        #layer_cols = pd.Index().difference(EXTRA_COLS)

        #if "read_id" in self.prms.tsv_cols:
        #    self.index_cols

        #self.extras = EXTRA_COLS.intersection(self.prms.tsv_cols)

        self.layers = INDEX_COLS+self.prms.tsv_cols+["seq.kmer"]
        self.columns = pd.MultiIndex.from_tuples(self.layers)
        self.conf.tracks.layers += self.layers
        self._header = True

        TrackIO._init_output(self, self.prms.buffered)

    def init_read_mode(self):
        raise RuntimeError("Reading from TSV not yet supported")

    def write_alignment(self, aln):
        df = aln.to_pandas(self.layers, INDEX_COLS).sort_index()
        
        for col in df.columns[df.columns.get_level_values(-1).str.endswith("kmer")]:
            kmers = df[col].dropna()
            df[col] = self.tracks.model.kmer_to_str(kmers)

        #for col in self.extras:
        #    df[col] = EXTRA_FNS[col](aln)

        df.reset_index(inplace=True, drop=self.prms.tsv_noref)

        if self._header:
            self.columns = df.columns
            labels = [
                ".".join([c for c in col if len(c) > 0]) 
                for col in self.columns]
            self._header = False
            if self.prms.buffered:
                self._set_output(tuple(labels))
            else:
                self._set_output("\t".join(labels) + "\n")

        df = df.reindex(columns=self.columns)#.loc[:,self.columns]

        out = df.to_csv(sep="\t", header=False, na_rep=self.prms.tsv_na, index=False, float_format="%.6g")

        self._set_output(out)

    def write_buffer(self, buf=None):
        if buf is None:
            buf = [self.out_buffer]

        for out in buf:
            if isinstance(out, tuple):
                if self._header:
                    self._header = False
                    self.output.write("\t".join(out) + "\n")
            else:
                self.output.write(out)

    def close(self):
        if not self.prms.buffered:
            self.output.close()

