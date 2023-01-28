import numpy as np
import pandas as pd
import sys
import pysam
import array
import os
import re
import time
from collections import defaultdict, deque
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from ...index import RefCoord
from ... import PoreModel
from ... import Config
from . import TrackIO
import _uncalled 

REF_TAG  = "ur"
SAMP_TAG = "us"
NORM_TAG = "un"
LENS_TAG = "ul"
CURS_TAG = "uc"
STDS_TAG = "ud"

class BAM(TrackIO):
    FORMAT = "bam"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

        if write:
            if not self.prms.buffered:
                if tracks.bam_in is None:
                    raise ValueError("No BAM template provided")
                self.header = tracks.bam_in.input.header.to_dict()
            else:
                self.header = self.prms.bam_header
            self.init_write_mode()
            self.input = None
        else:
            self.output = None
            self.init_read_mode()

    def init_read_mode(self):
        if self.conf.tracks.io.bam_in is None:
            raise ValueError("BAM output is only available if BAM input provided")
        self.input = pysam.AlignmentFile(self.filename, "rb")

        self.header = self.input.header.to_dict()

        conf = None
        if "CO" in self.header:
            for line in self.header["CO"]:
                if not line.startswith("UNC:"): continue
                toml = re.sub(r"(?<!\\);", "\n", line[4:]).replace("\\;",";")
                conf = Config(toml=toml)

        if conf is None:
            conf = self.conf
            #toml = self.conf.to_toml()

        name = os.path.basename(self.filename)
        self.in_id = self.init_track(name, name, conf)

        self.in_alns = None

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

        if self.prms.buffered:
            self.out_buffer = list()
            self.output = None
            return

        #Store config toml in single-line comment with newlines replaced by semicolons
        conf_line = self.conf.to_toml() \
                        .replace(";", "\\;") \
                        .replace("\n", ";")   
        #header = self.input.header.to_dict()
        if not "CO" in self.header:
            self.header["CO"] = list()
        self.header["CO"].append("UNC:" + conf_line)

        self.output = pysam.AlignmentFile(self.conf.tracks.io.bam_out, "wb", header=self.header)#template=self.input)

    def get_alns(self, read_id):
        self._init_alns()
        return self.in_alns[read_id]

    def get_aln(self, read_id, ref_name, ref_start):
        self._init_alns()
        for aln in self.in_alns[read_id]:
            if aln.reference_name == ref_name and aln.reference_start == ref_start:
                return aln
        return None

    INF_U16 = np.iinfo(np.uint16).max

    def _init_alns(self):
        if self.in_alns is None:
            self.in_alns = defaultdict(list)
            for aln in self.iter_sam():
                self.in_alns[aln.query_name].append(aln)
            self.input.reset()

    def write_layers(self, track, groups):
        aln = track.alignments.iloc[0]
        #sam = self.get_aln(aln["read_id"], aln["ref_name"], aln["ref_start"]-self.kmer_trim[0])
        if self.bam is None: 
            return
        
        #refs = track.layer_refs
        #if aln["fwd"]:
        refs = track.coords.refs[self.kmer_trim[0]:-self.kmer_trim[1]]
        #else:
        #    refs = track.coords.refs[self.kmer_trim[1]:-self.kmer_trim[0]]

        dtw = track.layers["dtw"].reset_index(level=1).reindex(refs)

        lens = dtw["length"]
        lens[dtw["start"].duplicated()] = 0 #skips

        cmin = 0
        cmax = dtw["current"].max()
        scale = (self.INF_U16 - 1) / (cmax - cmin)
        shift = -cmin

        dc = np.round((dtw["current"]+ shift) * scale)#.astype("uint16")
        ds = np.round((dtw["stdv"]+ shift) * scale)#.astype("uint16")

        dc = dc.fillna(self.INF_U16).to_numpy(np.uint16)
        ds = ds.fillna(self.INF_U16).to_numpy(np.uint16)
        lens = dtw["length"].reindex(refs).to_numpy(np.uint16, na_value=0)

        if dtw["start"].iloc[0] < dtw["start"].iloc[-1]:
            sample_start = dtw["start"].iloc[0]
            sample_end = dtw["start"].iloc[-1] + dtw["length"].iloc[-1]
        else:
            sample_start = dtw["start"].iloc[0] + dtw["length"].iloc[0]
            sample_end = dtw["start"].iloc[-1] 

        self.bam.set_tag(REF_TAG, array.array("I", (refs[0], refs[-1]+1)))
        self.bam.set_tag(SAMP_TAG, array.array("I", (sample_start, sample_end)))
        self.bam.set_tag(NORM_TAG, array.array("f", (scale,shift)))
        self.bam.set_tag(LENS_TAG, array.array("H", lens))
        self.bam.set_tag(CURS_TAG, array.array("H", dc))
        self.bam.set_tag(STDS_TAG, array.array("H", ds))

        if not self.bam.has_tag("f5"):
            self.bam.set_tag("f5", os.path.basename(self.prev_fast5[0]))

        if self.prms.buffered:
            self.out_buffer.append(self.bam.to_string())
        else:
            self.output.write(self.bam)

    def write_buffer(self, bams):
        header = pysam.AlignmentHeader.from_dict(self.header)
        for b in bams:
            bam = pysam.AlignedSegment.fromstring(b, header)
            self.output.write(bam)

    #TODO more query options
    def iter_sam(self, unmapped=False):
        #if self.conf.tracks.ref_bounds is not None
        if self.conf.tracks.ref_bounds is None:
            itr = self.input
        else:
            b = self.conf.tracks.ref_bounds
            itr = self.input.fetch(b.name, b.start, b.end)

        if unmapped:
            mapped = lambda a: True
        else:
            mapped = lambda a: not a.is_unmapped

        read_filter = self.tracks.read_index.read_filter
        if read_filter is not None:
            filt = lambda a: a.query_name in read_filter
        else:
            filt = lambda a: True

        valid = lambda a: mapped(a) and filt(a)
            
        if self.conf.tracks.max_reads is None:
            for a in itr:
                if valid(a): yield a
        else:
            n = 0
            for a in itr:
                if valid(a): 
                    yield a
                    n += 1
                if n == self.conf.tracks.max_reads:
                    break
            

    def _parse_sam(self, sam, coords=None):
        samp_bounds = sam.get_tag(SAMP_TAG)

        norm_scale, norm_shift = sam.get_tag(NORM_TAG)

        dac = np.array(sam.get_tag(CURS_TAG))
        current = (dac/norm_scale)-norm_shift
        current[dac == self.INF_U16] = np.nan

        das = np.array(sam.get_tag(STDS_TAG))
        stdv = (das/norm_scale)-norm_shift
        stdv[das == self.INF_U16] = np.nan

        
        length = sam.get_tag(LENS_TAG)

        fwd = int(not sam.is_reverse)
        ref_start, ref_end = sam.get_tag(REF_TAG)


        if coords is None:
            refs = pd.RangeIndex(ref_start, ref_end)
            seq_refs = RefCoord(sam.reference_name, ref_start, ref_end, fwd)
            coords = self.tracks.index.get_coord_space(seq_refs, self.conf.is_rna, load_kmers=True)
            start_shift = end_shift = 0
        else:
            start_clip = (coords.refs.min() - ref_start) if coords.refs.min() > ref_start else 0
            end_clip = (ref_end - coords.refs.max()) if coords.refs.max() < ref_end else 0

            new_end = len(current) - end_clip
            current = current[start_clip:new_end]
            stdv = stdv[start_clip:new_end]


            start_shift = np.sum(length[:start_clip])
            end_shift = np.sum(length[new_end:])
            length = length[start_clip:new_end]

            ref_start += start_clip
            ref_end -= end_clip
            refs = pd.RangeIndex(ref_start, ref_end)

        pacs = coords.ref_to_pac(refs)
            
        kmers = coords.ref_kmers.loc[fwd].loc[refs].to_numpy()

        dtw = pd.DataFrame({
            "track_id" : self.in_id, "fwd" : fwd, "pac" : pacs, "aln_id" : self.aln_id_in,
            "current" : current, "stdv" : stdv, "length" : length, "kmer" : kmers
        }).set_index(["track_id", "fwd", "pac","aln_id"])

        dtw.loc[dtw["current"].isnull(), "kmer"] = pd.NA
        dtw["kmer"] = dtw["kmer"].astype("Int32", copy=False)


        flip = samp_bounds[0] > samp_bounds[1]
        if flip:
            samp_end, samp_start = samp_bounds
            samp_start += end_shift
            samp_end -= start_shift
            dtw = dtw.iloc[::-1]
        else:
            samp_start, samp_end = samp_bounds
            samp_start += start_shift
            samp_end -= end_shift

        start = np.full(len(dtw), samp_start, dtype="int32")
        start[1:] += dtw["length"].cumsum().iloc[:-1].to_numpy("int32")
        dtw["start"] = start

        if flip:
            dtw = dtw[::-1]

        layers = pd.concat({"dtw" : dtw}, names=("group","layer"), axis=1)#.sort_index()

        layers.loc[layers["dtw", "length"] == 0, ("dtw", "length")] = pd.NA
        layers = layers[~(layers["dtw","current"].isnull() & layers["dtw","length"].isnull())]

        layers["dtw", "length"] = layers["dtw", "length"].fillna(method="pad").astype("Int32")

        #Note should use below, but too buggy: https://github.com/pandas-dev/pandas/issues/45725
        #layers["dtw", "length"] = layers["dtw", "length"].replace(0, pd.NA)
        #layers["dtw", "length"] = layers["dtw", "length"].fillna(method="pad").astype("int32")

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
        
        new_layers = list()
        new_alns = list()

        layer_df = pd.DataFrame()
        aln_df = pd.DataFrame()

        prev_ref = None
        prev_start = 0
        #########prev_aln = 0

        t = time.time()
        for sam in itr:
            next_alns,next_layers = self._parse_sam(sam)

            next_aln = next_alns.iloc[0]
            next_ref = next_aln["ref_name"]
            next_start = next_layers.index.get_level_values("pac").min()

            #we're losing alignments somewhere. getting out of sync
            #maybe right below
            #but also, not sure if I'm syncing inputs in "Tracks"
            #probably need to always advance the minimum ref pos/ID

            ret_layer_count = 0

            if prev_ref is not None and (next_ref != prev_ref or next_start - prev_start > 0) and len(new_alns) > self.prms.aln_chunksize:
                #pac = self.tracks.index.ref_to_pac(a["ref_name"], prev_start)

                aln_df = pd.concat([aln_df] + new_alns)#.sort_index()
                layer_df = pd.concat([layer_df] + new_layers).sort_index()
                new_alns = list()
                new_layers = list()
                
                ret_layers = layer_df.loc[slice(None),slice(None),prev_start:next_start-1]
                ret_alns = aln_df.loc[ret_layers.index.droplevel(["fwd", "pac"]).unique()]

                layer_df = layer_df.loc[slice(None),slice(None),next_start:]
                aln_df = aln_df.loc[layer_df.index.droplevel(["fwd", "pac"]).unique()]
                #ret_layer_count = 0

                t = time.time()

                if len(ret_alns) > 0:
                    yield (ret_alns, ret_layers)

                prev_ref = next_ref
                prev_start = next_start

            elif prev_ref is None:
                prev_ref = next_ref
                prev_start = next_start


            new_alns.append(next_alns)
            new_layers.append(next_layers)

        aln_df = pd.concat([aln_df] + new_alns).sort_index()
        layer_df = pd.concat([layer_df] + new_layers).sort_index()
        
        ret_layers = layer_df.loc[slice(None),slice(None),prev_start:]
        ret_alns = aln_df.loc[ret_layers.index.droplevel(["fwd", "pac"]).unique()]

        if len(ret_alns) > 0:
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
            ref_start, ref_end = sam.get_tag(REF_TAG)
            ref_bounds = RefCoord(sam.reference_name, ref_start, ref_end, not sam.is_reverse)
            coords = ref_index.get_coord_space(ref_bounds, is_rna=self.conf.is_rna, load_kmers=True, kmer_trim=False)

            yield self._parse_sam(sam, coords)

            aln_id += 1


    def close(self):
        if self.input is not None:
            self.input.close()
        if self.output is not None:
            self.output.close()

