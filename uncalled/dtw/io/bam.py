import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from . import TrackIO
import _uncalled 
import pysam
import array

class BAM(TrackIO):
    FORMAT = "bam"

    def __init__(self, conf, mode):
        self.in_fname = conf.tracks.io.bam_out
        TrackIO.__init__(self, conf.tracks.io.bam_in, conf, mode)

        self.init_read_mode()
        self.init_write_mode()

    def init_read_mode(self):
        if self.conf.tracks.io.bam_in is None:
            raise ValueError("BAM output is only available if BAM input provided")
        self.input = pysam.AlignmentFile(self.conf.tracks.io.bam_in, "rb")

        #TODO temporary until I replace PAF file with BAM input
        if self.output is not None:
            self.in_alns = {
                (aln.query_name, aln.reference_name, aln.reference_start+2) : aln
                for aln in self.input if not aln.is_unmapped
            }

    def init_write_mode(self):
        if self.conf.tracks.io.bam_out is None:
            return
        TrackIO.init_write_mode(self)
        self.output = pysam.AlignmentFile(self.conf.tracks.io.bam_out, "wb", template=self.input)

    def get_header(self, track):
        return {"HD" : {"VN" : "1.0"},
                "SQ" : [
                    {"SN" : name, "LN" : length}
                    for name, length in track.index.get_seqs()
                ]}

    INF_U16 = np.iinfo(np.uint16).max

    def write_layers(self, track, groups):
        aln = track.alignments.iloc[0]
        k = (aln["read_id"], aln["ref_name"], aln["ref_start"])
        if not k in self.in_alns: return
        sam = self.in_alns[k]

        dtw = track.layers["dtw"].dropna().reset_index(level=1)

        lens = dtw["length"]

        cmin = dtw["current"].min()
        cmax = dtw["current"].max()
        scale = (self.INF_U16 - 1) / (cmax - cmin)
        shift = -cmin

        dc = np.round((dtw["current"]+ shift) * scale)#.astype("uint16")

        refs = track.coords.refs[2:-2]
        dc = dc.reindex(refs, fill_value=self.INF_U16).to_numpy(np.uint16)
        lens = dtw["length"].reindex(refs, fill_value=self.INF_U16).to_numpy(np.uint16)

        #sam.set_tag("sc", dc.tobytes().hex(), "H") #hex format
        #sam.set_tag("sl", lens.tobytes().hex(), "H") #hex formt
        #cur = dtw["current"].to_numpy().astype(float)
        #sam.set_tag("sc", array.array("f", cur)) #float current
        #sam.set_tag("sl", ",".join(map(str,lens))) #str lens

        if dtw["start"].iloc[0] < dtw["start"].iloc[-1]:
            sample_start = dtw["start"].iloc[0]
            sample_end = dtw["start"].iloc[-1] + dtw["length"].iloc[-1]
        else:
            sample_start = dtw["start"].iloc[0] + dtw["length"].iloc[0]
            sample_end = dtw["start"].iloc[-1]

        sam.set_tag("sr", array.array("I", (refs[0], refs[-1])))
        sam.set_tag("ss", array.array("I", (sample_start, sample_end)))
        sam.set_tag("sn", array.array("f", (scale,shift)))
        sam.set_tag("sl", array.array("H", lens))
        sam.set_tag("sc", array.array("H", dc))

        self.output.write(sam)

        #df = track.layers_desc_index#.reset_index()
        #df = df[self.columns.intersection(df.columns)].dropna(how="all", axis=0).reset_index()

        #if self._header:
        #    self.columns = df.columns
        #    labels = [
        #        ".".join([c for c in col if len(c) > 0]) 
        #        for col in self.columns]
        #    self.out.write("\t".join(labels) + "\n")
        #    self._header = False

        #df = df.loc[:,self.columns]

        #self.out.write(df.to_csv(sep="\t", header=False, na_rep=self.prms.bam_na, index=False))

    def close(self):
        self.input.close()
        if self.output is not None:
            self.output.close()

