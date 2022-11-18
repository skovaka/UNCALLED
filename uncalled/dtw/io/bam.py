import numpy as np
import pandas as pd
import sys
import pysam
import array
import os
import re
from collections import defaultdict, deque
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from ...index import RefCoord
from ... import PoreModel
from ... import Config
from . import TrackIO
import _uncalled 

class BAM(TrackIO):
    FORMAT = "bam"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

        self.init_read_mode()
        self.init_write_mode()

    def init_read_mode(self):
        if self.conf.tracks.io.bam_in is None:
            raise ValueError("BAM output is only available if BAM input provided")
        self.input = pysam.AlignmentFile(self.filename, "rb")

        conf = None
        if "CO" in self.input.header:
            for line in self.input.header["CO"]:
                if not line.startswith("UNC:"): continue
                toml = re.sub(r"(?<!\\);", "\n", line[4:]).replace("\\;",";")
                conf = Config(toml=toml)

        if conf is None:
            conf = self.conf
            #toml = self.conf.to_toml()

        name = os.path.basename(self.filename)
        self.in_id = self.init_track(name, name, conf)

        #TODO really not good, should stream by default
        #if self.conf.tracks.io.bam_out is not None:
        self.in_alns = defaultdict(list)
        for aln in self.iter_sam():
            self.in_alns[aln.query_name].append(aln)
        self.input.reset()

        self.aln_id_in = 1

    def init_write_mode(self):
        if self.conf.tracks.io.bam_out is None:
            self.output = None
            return
        TrackIO.init_write_mode(self)

        if len(self.conf.fast5_reader.fast5_index) == 0:
            sys.stderr.write("Warning: no fast5 index specified\n")

        #TODO load from AlnTrack instance (initialized by Tracks)
        self.model = PoreModel(self.conf.pore_model)
        self.kmer_trim = self.model.kmer_trim

        #Store config toml in single-line comment with newlines replaced by semicolons
        conf_line = self.conf.to_toml() \
                        .replace(";", "\\;") \
                        .replace("\n", ";")   
        header = self.input.header.to_dict()
        if not "CO" in header:
            header["CO"] = list()
        header["CO"].append("UNC:" + conf_line)

        self.output = pysam.AlignmentFile(self.conf.tracks.io.bam_out, "wb", header=header)#template=self.input)

    def get_alns(self, read_id):
        return self.in_alns[read_id]

    def get_aln(self, read_id, ref_name, ref_start):
        for aln in self.in_alns[read_id]:
            if aln.reference_name == ref_name and aln.reference_start == ref_start:
                return aln
        return None

    INF_U16 = np.iinfo(np.uint16).max

    def write_layers(self, track, groups):
        aln = track.alignments.iloc[0]
        sam = self.get_aln(aln["read_id"], aln["ref_name"], aln["ref_start"]-self.kmer_trim[0])
        if sam is None: 
            return

        refs = track.coords.refs[self.kmer_trim[0]:-self.kmer_trim[1]]

        dtw = track.layers["dtw"].reset_index(level=1).reindex(refs)

        lens = dtw["length"]
        lens[dtw["start"].duplicated()] = 0 #skips

        cmin = dtw["current"].min()
        cmax = dtw["current"].max()
        scale = (self.INF_U16 - 1) / (cmax - cmin)
        shift = -cmin

        dc = np.round((dtw["current"]+ shift) * scale)#.astype("uint16")

        dc = dc.fillna(self.INF_U16).to_numpy(np.uint16)
        lens = dtw["length"].reindex(refs).to_numpy(np.uint16, na_value=0)

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

    #TODO more query options
    def iter_sam(self, unmapped=False):
        #if self.conf.tracks.ref_bounds is not None
        for sam in self.input:
            if unmapped or not sam.is_unmapped:
                yield sam

    def _parse_sam(self, sam, coords=None):
        samp_bounds = sam.get_tag("ss")

        norm_scale, norm_shift = sam.get_tag("sn")
        dac = np.array(sam.get_tag("sc"))
        current = (dac/norm_scale)-norm_shift
        current[dac == self.INF_U16] = np.nan
        
        length = sam.get_tag("sl")


        fwd = int(not sam.is_reverse)
        ref_start, ref_end = sam.get_tag("sr")
        refs = pd.RangeIndex(ref_start, ref_end)

        if coords is None:
            seq_refs = RefCoord(sam.reference_name, ref_start, ref_end, fwd)
            coords = self.tracks.index.get_coord_space(seq_refs, self.conf.is_rna, load_kmers=True)
            
        pacs = coords.ref_to_pac(refs)
        kmers = coords.ref_kmers.loc[fwd].loc[refs].to_numpy()

        dtw = pd.DataFrame(
            {"current" : current, "length" : length},
            index=pd.MultiIndex.from_product([[self.in_id], [fwd], pacs, [self.aln_id_in]],
                                             names=("track_id", "fwd","pac","aln_id"))
        )

        if samp_bounds[0] < samp_bounds[1]:
            samp_start, samp_end = samp_bounds
        else:
            samp_end, samp_start = samp_bounds
            dtw = dtw.iloc[::-1]

        start = np.full(len(dtw), samp_start, dtype="int32")
        start[1:] += dtw["length"].cumsum().iloc[:-1].to_numpy("int32")
        dtw["start"] = start

        layers = pd.concat({"dtw" : dtw}, names=("group","layer"), axis=1).sort_index()

        layers["dtw", "length"] = layers["dtw", "length"].replace(0, pd.NA)#.astype("int32")
        layers["dtw", "length"] = layers["dtw", "length"].fillna(method="pad").astype("int32")

        aln = pd.DataFrame({
                "id" : [self.aln_id_in],
                "track_id" : [self.in_id],
                "read_id" : [sam.query_name],
                "ref_name" : [coords.ref_name],
                "ref_start" : [ref_start],
                "ref_end" : [ref_end],
                "fwd" :     [fwd],
                "samp_start" : [samp_start],
                "samp_end" : [samp_end],
                "tags" : [""]}).set_index(["track_id", "id"])

        self.aln_id_in += 1

        return aln, layers

    def query(self, layers, track_id, coords, aln_id=None, read_id=None, fwd=None, order=["read_id", "pac"], full_overlap=False):
        itr = self.input.fetch(coords.ref_name, coords.refs.min(), coords.refs.max())

        alignments = list()
        layers = list()

        if read_id is not None and len(read_id) > 0:
            read_ids = set(read_id)
        else:
            read_ids = None

        for sam in itr:
            if read_ids is not None and sam.query_name is not read_ids:
                continue
            a,l = self._parse_sam(sam, coords)
            alignments.append(a)
            layers.append(l)

        return pd.concat(alignments), pd.concat(layers)
        #self.fill_tracks(coords, pd.concat(alignments), pd.concat(layers))

    def iter_refs(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, chunksize=None, full_overlap=False):

        if coords is not None:
            itr = self.input.fetch(coords.ref_name, coords.refs.min(), coords.refs.max())
        else:
            itr = self.input.fetch()

        if fwd is not None:
            strands = [int(fwd)]
        else:
            strands = [1, 0]
        
        layers = deque()
        alns = deque()

        prev_ref = None
        prev_start = 0
        prev_aln = 0

        for sam in itr:
            aln,layer = self._parse_sam(sam)

            a = aln.iloc[0]
            pac_st = layer.index.get_level_values("pac").min()

            if prev_ref is not None and (a["ref_name"] != prev_ref or a["ref_start"] > prev_start) and len(alns) > self.prms.aln_chunksize:
                pac = self.tracks.index.ref_to_pac(a["ref_name"], prev_start)
                ret_layers = pd.concat([l[l.index.get_level_values("pac").isin(pd.RangeIndex(prev_start, pac_st))] for l in layers])
                ret_alns = pd.concat([alns[prev_aln:]])

                yield (ret_alns, ret_layers)

                while len(alns) > 0 and (alns[0].iloc[0]["ref_name"] != a["ref_name"] or layers[0].index.get_level_values("pac").max() < pac_st):
                    alns.popleft()
                    layers.popleft()

                prev_ref = a["ref_name"]
                prev_start = a["ref_start"]
                prev_aln = len(alns)

            alns.append(aln)
            layers.append(layer)

        pac = self.tracks.index.ref_to_pac(a["ref_name"], prev_start)
        ret_layers = pd.concat([l[l.index.get_level_values("pac") >= prev_start] for l in layers])
        #aln_ids = ret_layers.index.get_level_values("aln_id").unique()
        aln_ids = ret_layers.index.droplevel(["fwd","pac"]).unique()
        ret_alns = pd.concat([a.loc[aln_ids.intersection(a.index)] for a in alns])

        yield (ret_alns, ret_layers)

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
        aln_id = 0

        if read_id is not None and len(read_id) > 0:
            read_ids = set(read_id)
        else:
            read_ids = None

        for sam in self.iter_sam():
            if read_ids is not None and sam.query_name not in read_ids:
                continue
            ref_start, ref_end = sam.get_tag("sr")
            ref_bounds = RefCoord(sam.reference_name, ref_start, ref_end, not sam.is_reverse)
            coords = ref_index.get_coord_space(ref_bounds, is_rna=self.conf.is_rna, load_kmers=True, kmer_trim=False)

            ##refs = pd.RangeIndex(*sam.get_tag("sr"), name="ref")
            #samp_bounds = sam.get_tag("ss")

            #norm_scale, norm_shift = sam.get_tag("sn")
            #dac = np.array(sam.get_tag("sc"))
            #current = (dac/norm_scale)-norm_shift
            #current[dac == self.INF_U16] = np.nan
            #
            #length = sam.get_tag("sl")

            #dtw = pd.DataFrame(
            #    {"current" : current, "length" : length, "kmer" : coords.ref_kmers.to_numpy()},
            #    index=pd.MultiIndex.from_product([[int(coords.fwd)], coords.pacs, [aln_id]], 
            #                                     names=("fwd","pac","aln_id"))
            #)

            #
            #if samp_bounds[0] < samp_bounds[1]:
            #    samp_start, samp_end = samp_bounds
            #else:
            #    samp_end, samp_start = samp_bounds
            #    dtw = dtw.iloc[::-1]

            #start = np.full(len(dtw), samp_start, dtype="int32")
            #start[1:] += dtw["length"].cumsum().iloc[:-1].to_numpy("int32")
            #dtw["start"] = start

            #layers = pd.concat({"dtw" : dtw}, names=("group","layer"), axis=1).sort_index()

            #alns = pd.DataFrame({
            #        "id" : [aln_id],
            #        "track_id" : [1],
            #        "read_id" : [sam.query_name],
            #        "ref_name" : [coords.ref_name],
            #        "ref_start" : [coords.refs.start],
            #        "ref_end" : [coords.refs.stop],
            #        "fwd" :     [int(coords.fwd)],
            #        "samp_start" : [samp_start],
            #        "samp_end" : [samp_end],
            #        "tags" : [""]}).set_index("id")

            yield self._parse_sam(sam, coords)

            aln_id += 1


    def close(self):
        self.input.close()
        if self.output is not None:
            self.output.close()

