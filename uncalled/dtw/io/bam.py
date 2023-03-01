import numpy as np
import pandas as pd
import sys
import pysam
import array
import os
import re
import time
from collections import defaultdict, deque
from ..bcaln import sam_to_ref_moves, INT32_NA
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from ...index import RefCoord
from ... import PoreModel
from ... import Config
from ...dfs import Alignment, AlnDF
from . import TrackIO
from _uncalled import IntervalIndexI64


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

        self._init_tags(conf.tracks.io.bam_tags)

        name = os.path.basename(self.filename)
        self.in_id = self.init_track(name, name, conf)

        self.in_alns = None

        self.aln_id_in = 1

    def reset(self):
        self.input.reset()
        self.aln_id_in = 1

    def _init_tags(self, tags):
        self.tags = dict()
        for t in tags:
            label, tag = t.split(":")
            self.tags[label] = tag

    def init_write_mode(self):
        if self.conf.tracks.io.bam_out is None:
            self.output = None
            return

        self._init_tags(self.prms.bam_tags)

        TrackIO.init_write_mode(self)

        #TODO load from AlnTrack instance (initialized by Tracks)
        self.model = PoreModel(self.conf.pore_model)
        #self.kmer_trim = self.model.kmer_trim

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

    def write_alignment(self, aln):
        if self.bam is None: 
            return

        refs = list()
        for i in range(aln.seq.coords.interval_count()):
            c = aln.seq.coords.get_interval(i)
            refs += [c.start, c.end]

        cmin = 0
        cmax = np.max(aln.dtw.current)
        scale = (self.INF_U16 - 1) / (cmax - cmin)
        shift = -cmin

        dc = np.round((aln.dtw.current.to_numpy() + shift) * scale).astype(np.uint16)
        ds = np.round((aln.dtw.current_sd.to_numpy() + shift) * scale).astype(np.uint16)

        lens = aln.dtw.samples.lengths_dedup.to_numpy()

        self.bam.set_tag(REF_TAG, array.array("I", refs))
        self.bam.set_tag(SAMP_TAG, array.array("I", (aln.dtw.samples.start, aln.dtw.samples.end)))
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

    def write_layers(self, track, groups):
        aln = track.alignments.iloc[0]
        if self.bam is None: 
            return
        
        refs = track.coords.refs
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

    def iter_str_chunks(self):
        read_ids = set()
        bams = list()
        for bam in self.iter_sam():
            read_ids.add(bam.query_name)
            bams.append(bam.to_string())

            if len(bams) == self.prms.bam_chunksize:
                yield(read_ids, bams)
                read_ids = set()
                bams = list()

        if len(bams) > 0:
            yield(read_ids, bams)

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

    def iter_moves(self):
        for sam in self.iter_sam():
            if sam.query_name in self.tracks.read_index:
                read = self.tracks.read_index[sam.query_name]

            moves = sam_to_ref_moves(self.conf, self.tracks.index, read, sam, self.tracks.read_index.default_model)
            if moves is None:
                continue
            self.bam = sam

            seq = self.tracks.get_seq(moves.index)
            aln = Alignment(read, seq, sam)
            aln.set_moves(moves)
            
            yield sam, aln

    def _parse_sam(self, sam, coords=None):
        samp_bounds = sam.get_tag(SAMP_TAG)

        norm_scale, norm_shift = sam.get_tag(NORM_TAG)

        dac = np.array(sam.get_tag(CURS_TAG))
        current = (dac/norm_scale)-norm_shift
        current[dac == self.INF_U16] = np.nan

        if sam.has_tag(STDS_TAG):
            das = np.array(sam.get_tag(STDS_TAG))
            stdv = (das/norm_scale)-norm_shift
            stdv[das == self.INF_U16] = np.nan
        else:
            stdv = np.full(len(current), np.nan)

        
        length = sam.get_tag(LENS_TAG)

        fwd = int(not sam.is_reverse)
        #ref_start, ref_end = sam.get_tag(REF_TAG)
        refs = sam.get_tag(REF_TAG)

        icoords = IntervalIndexI64([(refs[i], refs[i+1]) for i in range(0, len(refs), 2)])

        samp_start, samp_end = samp_bounds

        #TODO build whole alignment, slice after
        if coords is not None:
            clip = IntervalIndexI64([(coords.mrefs[fwd].min(), coords.mrefs[fwd].max()+1)])

            st = icoords.get_index(clip.start) if icoords.start < clip.start else 0
            en = icoords.get_index(clip.end-1)+1  if icoords.end > clip.end else icoords.length

            icoords = icoords.islice(st,en)

            current = current[st:en]
            stdv = stdv[st:en]

            samp_start += np.sum(length[:st])
            samp_end -= np.sum(length[en:])
            length = length[st:en]

        pacs = self.tracks.index.mref_to_pac(icoords.expand())

        seq = self.tracks.get_seq(icoords)
        read = self.tracks.read_index.get(sam.query_name, None)
        aln = Alignment(sam.query_name if read is None else read, seq, sam)
            
        kmers = seq.kmer #coords.ref_kmers.loc[fwd].loc[refs].to_numpy()

        start = np.full(icoords.length, samp_start, dtype="int32")
        start[1:] += np.cumsum(length)[:-1].astype("int32")
        #.to_numpy("int32")

        dtw = AlnDF(icoords, start, length, current, stdv)
        aln.set_dtw(dtw)

        moves = sam_to_ref_moves(self.conf, self.tracks.index, read, sam, self.tracks.read_index.default_model)
        if moves is not None:
            aln.set_moves(moves)

        dtw = pd.DataFrame({
            "track_id" : self.in_id, "fwd" : fwd, "pac" : pacs, "aln_id" : self.aln_id_in,
            "current" : current, "stdv" : stdv, "start" : start, "length" : length, "kmer" : seq.kmer
        }).set_index(["track_id", "fwd", "pac","aln_id"])

        dtw.loc[dtw["current"].isnull(), "kmer"] = pd.NA
        dtw["kmer"] = dtw["kmer"].astype("Int32", copy=False)

        layers = pd.concat({"dtw" : dtw}, names=("group","layer"), axis=1)#.sort_index()

        layers.loc[layers["dtw", "length"] == 0, ("dtw", "length")] = pd.NA
        layers = layers[~(layers["dtw","current"].isnull() & layers["dtw","length"].isnull())]

        layers["dtw", "length"] = layers["dtw", "length"].fillna(method="pad").astype("Int32")

        aln_df = pd.DataFrame({
                "id" : [self.aln_id_in],
                "track_id" : [self.in_id],
                "read_id" : [sam.query_name],
                "ref_name" : [sam.reference_name],
                "ref_start" : [sam.reference_start],
                "ref_end" : [sam.reference_end],
                "fwd" :     [fwd],
                "samp_start" : [samp_start],
                "samp_end" : [samp_end],
                "tags" : [""]}).set_index(["track_id", "id"])

        if moves is not None:
            pacs = self.tracks.index.mref_to_pac(moves.index.expand())
            df = pd.DataFrame({
                #"mref"   : (moves.index.expand()),
                "track_id" : self.in_id, "fwd" : fwd, "pac" : pacs, "aln_id" : self.aln_id_in,
                "start"  : (moves.samples.starts),
                "length" : (moves.samples.lengths)
            }).set_index(["track_id", "fwd", "pac","aln_id"])
            #}).set_index("mref")
            isna = df["start"] == INT32_NA
            df[isna] = pd.NA
            df = pd.concat({"bcaln" : df.fillna(method="backfill")}, names=("group","layer"),  axis=1)

            layers = pd.concat([layers, df.reindex(dtw.index)], axis=1)

        aln.calc_mvcmp()
        pacs = self.tracks.index.mref_to_pac(aln._mvcmp.index.expand())
        cmp_df = pd.DataFrame({
                "track_id" : self.in_id, "fwd" : fwd, "pac" : pacs, "aln_id" : self.aln_id_in,
                "aln_b" : self.aln_id_in, "mean_ref_dist" : aln._mvcmp.dist, "jaccard" : aln._mvcmp.jaccard,
        }).set_index(["track_id", "fwd", "pac","aln_id"])
        #print(cmp_df)
        #print(layers)

        layers = pd.concat([layers, pd.concat({"bc_cmp" :cmp_df}, axis=1)], axis=1)

        self.aln_id_in += 1

        #print("dtw")
        #print(aln.to_pandas())
        #print(layers)

        return aln_df, layers

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

    def iter_alns(self, layers=None, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
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

