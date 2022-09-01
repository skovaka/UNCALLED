import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from ...index import RefCoord
from . import TrackIO
import _uncalled 
import pysam
import array
import os

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
        filename = self.conf.tracks.io.bam_in
        self.input = pysam.AlignmentFile(filename, "rb")

        name = os.path.basename(filename)
        self.init_track(1, name, name, self.conf.to_toml())

        #TODO temporary until I replace PAF file with BAM input
        if self.conf.tracks.io.bam_out is not None:
            self.in_alns = {
                (aln.query_name, aln.reference_name, aln.reference_start+2) : aln
                for aln in self.input if not aln.is_unmapped
            }

    def init_write_mode(self):
        if self.conf.tracks.io.bam_out is None:
            self.output = None
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

        refs = track.coords.refs[2:-2]

        dtw = track.layers["dtw"].reset_index(level=1).reindex(refs)

        lens = dtw["length"]
        lens[dtw["start"].duplicated()] = 0 #skips

        cmin = dtw["current"].min()
        cmax = dtw["current"].max()
        scale = (self.INF_U16 - 1) / (cmax - cmin)
        shift = -cmin

        dc = np.round((dtw["current"]+ shift) * scale)#.astype("uint16")

        dc = dc.fillna(self.INF_U16).to_numpy(np.uint16)
        lens = dtw["length"].reindex(refs, fill_value=0).to_numpy(np.uint16)

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

        sam.set_tag("sr", array.array("I", (refs[0], refs[-1]+1)))
        sam.set_tag("ss", array.array("I", (sample_start, sample_end)))
        sam.set_tag("sn", array.array("f", (scale,shift)))
        sam.set_tag("sl", array.array("H", lens))
        sam.set_tag("sc", array.array("H", dc))

        self.output.write(sam)

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
        
        aln_id = 0

        for sam in self.input:
            ref_start, ref_end = sam.get_tag("sr")
            ref_bounds = RefCoord(sam.reference_name, ref_start, ref_end, not sam.is_reverse)
            coords = ref_index.get_coord_space(ref_bounds, is_rna=self.conf.is_rna, load_kmers=True, kmer_trim=False)

            #refs = pd.RangeIndex(*sam.get_tag("sr"), name="ref")
            samp_bounds = sam.get_tag("ss")

            norm_scale, norm_shift = sam.get_tag("sn")
            dac = np.array(sam.get_tag("sc"))
            current = (dac/norm_scale)-norm_shift
            current[dac == self.INF_U16] = np.nan
            
            length = sam.get_tag("sl")

            dtw = pd.DataFrame(
                {"current" : current, "length" : length, "kmer" : coords.ref_kmers.to_numpy()},
                index=pd.MultiIndex.from_product([[int(coords.fwd)], coords.pacs, [aln_id]], 
                                                 names=("fwd","pac","aln_id"))
            )

            
            if samp_bounds[0] < samp_bounds[1]:
                samp_start, samp_end = samp_bounds
            else:
                samp_end, samp_start = samp_bounds
                dtw = dtw.iloc[::-1]

            start = np.full(len(dtw), samp_start, dtype="int32")
            start[1:] += dtw["length"].cumsum().iloc[:-1].to_numpy("int32")
            dtw["start"] = start#[::-1]

            #start = np.full(len(dtw), samp_start, dtype="int32")
            #start[1:] += dtw["length"].cumsum().iloc[:-1].to_numpy("int32")


            layers = pd.concat({"dtw" : dtw}, names=("group","layer"), axis=1).sort_index()

            #def make_groups(df):
            #    grouped = dict()
            #    for group, layers in group_layers.items():
            #        gdf = df[layers]
            #        grouped[group] = df[layers].rename(columns=renames)
            #    df = pd.concat(grouped, names=("group", "layer"), axis=1)
            #    df.index.names = ("fwd", "pac", "aln_id")

            alns = pd.DataFrame({
                    "id" : [aln_id],
                    "track_id" : [1],
                    "read_id" : [sam.query_name],
                    "ref_name" : [coords.ref_name],
                    "ref_start" : [coords.refs.start],
                    "ref_end" : [coords.refs.stop],
                    "fwd" :     [int(coords.fwd)],
                    "samp_start" : [samp_start],
                    "samp_end" : [samp_end],
                    "tags" : [""]}).set_index("id")

            yield alns, layers

            aln_id += 1


    def close(self):
        self.input.close()
        if self.output is not None:
            self.output.close()

