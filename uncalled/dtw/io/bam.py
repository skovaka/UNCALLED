import numpy as np
import pandas as pd
import sys
import pysam
import array
import os
import re
import time
from collections import defaultdict, deque
from ..moves import sam_to_ref_moves, INT32_NA
from ..aln_track import AlnTrack
from ..layers import LAYER_META, parse_layers
from ...index import RefCoord
from ... import PoreModel
from ... import Config
from ...aln import Alignment, AlnDF
from . import TrackIO
from _uncalled import IntervalIndexI64, ValArrayI16


REF_TAG  = "ur"
SAMP_TAG = "us"
NORM_TAG = "un"
LENS_TAG = "ul"
CURS_TAG = "uc"
STDS_TAG = "ud"

REQ_ALN_TAGS = [REF_TAG, SAMP_TAG, LENS_TAG, CURS_TAG]

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

        self.read_id_index = None

        self.in_alns = None

    def reset(self):
        self.input.reset()

    def iter_read(self, read_id):
        if self.read_id_index is None:
            sys.stderr.write("Indexing BAM file by-read...\n")
            self.read_id_index = pysam.IndexedReads(self.input)
            self.read_id_index.build()
        ret = self.read_id_index.find(read_id)
        return ret


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
        self.model = PoreModel(params=self.conf.pore_model)
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

    INF_U16 = np.iinfo(np.int16).max
    NA_I16 = ValArrayI16.NA

    def _init_alns(self):
        if self.in_alns is None:
            self.in_alns = defaultdict(list)
            for aln in self.iter_sam():
                self.in_alns[aln.query_name].append(aln)
            self.input.reset()

    def write_alignment(self, aln):
        self.bam = aln.sam

        refs = list()
        for i in range(aln.seq.coords.interval_count()):
            c = aln.seq.coords.get_interval(i)
            refs += [c.start, c.end]

        cmin = np.min(aln.dtw.current)
        cmax = np.max(aln.dtw.current)
        #scale = (self.INF_U16 - 1) / (cmax - cmin)
        #shift = -cmin
        scale = self.model.current.INORM_SCALE
        shift = 0

        dc = np.round((aln.dtw.current.to_numpy() + shift) * scale)#.astype(np.int16)
        ds = np.round((aln.dtw.current_sd.to_numpy() + shift) * scale)#.astype(np.int16)

        na = np.isnan(dc)
        dc[na] = self.NA_I16
        ds[na] = self.NA_I16

        lens = aln.dtw.samples.lengths_dedup.to_numpy().astype("int16")

        self.bam.set_tag(REF_TAG, array.array("i", refs))
        self.bam.set_tag(SAMP_TAG, array.array("i", (aln.dtw.samples.start, aln.dtw.samples.end)))
        self.bam.set_tag(NORM_TAG, array.array("f", (scale,shift)))
        self.bam.set_tag(LENS_TAG, array.array("h", lens.astype(np.int16)))
        self.bam.set_tag(CURS_TAG, array.array("h", dc.astype(np.int16)))
        self.bam.set_tag(STDS_TAG, array.array("h", ds.astype(np.int16)))

        #if not self.bam.has_tag("f5"):
        #    self.bam.set_tag("f5", os.path.basename(self.prev_fast5[0]))

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

    def iter_alns(self, layers=None, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
        for sam in self.iter_sam():
            if sam.query_name in self.tracks.read_index:
                read = self.tracks.read_index[sam.query_name]

            moves = sam_to_ref_moves(self.conf, self.tracks.index, read, sam, self.tracks.read_index.default_model)
            if moves is None:
                continue
            self.bam = sam

            aln = self.tracks.init_alignment(self.next_aln_id(), read, sam.reference_id, moves.index, sam)
            aln.set_moves(moves)

            yield sam, aln

    def sam_to_aln(self, sam, load_moves=True):
        has_dtw = True
        for tag in REQ_ALN_TAGS:
            if not sam.has_tag(tag):
                has_dtw = False
                break

        aln = None

        try:
            read = self.tracks.read_index.get(sam.query_name, None)
        except:
            sys.stderr.write(f"Warning: failed to open read {sam.query_name}\n")
            read = None

        if not has_dtw and read is None:
            return None

        if has_dtw:
            samp_bounds = sam.get_tag(SAMP_TAG)

            norm_scale, norm_shift = sam.get_tag(NORM_TAG)

            dac = np.array(sam.get_tag(CURS_TAG))
            cna = dac == self.NA_I16
            current = (dac/norm_scale)-norm_shift
            current[cna] = np.nan

            if sam.has_tag(STDS_TAG):
                das = np.array(sam.get_tag(STDS_TAG))
                stdv = (das/norm_scale)-norm_shift
                stdv[das == self.NA_I16] = np.nan
            else:
                stdv = np.full(len(current), np.nan)
            
            length = np.array(sam.get_tag(LENS_TAG))
            #ref_start, ref_end = sam.get_tag(REF_TAG)
            refs = sam.get_tag(REF_TAG)

            coords = IntervalIndexI64([(refs[i], refs[i+1]) for i in range(0, len(refs), 2)])

            samp_start, samp_end = samp_bounds

            aln = self.tracks.init_alignment(self.next_aln_id(), read, sam.reference_id, coords, sam)

            start = np.full(coords.length, samp_start, dtype="int32")
            start[1:] += np.cumsum(length.astype("int32"))[:-1]

            lna = (length == 0) & cna
            length[cna] = -1

            dtw = AlnDF(aln.seq, start, length, current, stdv)
            aln.set_dtw(dtw)



        moves = sam_to_ref_moves(self.conf, self.tracks.index, read, sam, self.tracks.read_index.default_model)
        has_moves = moves is not None
        if has_moves and load_moves:
            if aln is None:
                aln = self.tracks.init_alignment(self.next_aln_id(), read, sam.reference_id, moves.index, sam)
            else:
                i = max(0, moves.index.start - aln.seq.coords.start)
                j = min(len(moves), len(moves) + moves.index.end -  aln.seq.coords.end)
                moves = moves.slice(i,j)
                
            aln.set_moves(moves)

        if has_moves and has_dtw:
            aln.calc_mvcmp()
        elif aln is None:
            fwd = int(not sam.is_reverse)
            if fwd == self.conf.is_rna:
                coords = IntervalIndexI64([(-sam.reference_end, -sam.reference_start)])
            else:
                coords = IntervalIndexI64([(sam.reference_start, sam.reference_end)])
            aln = self.tracks.init_alignment(self.next_aln_id(), read, sam.reference_id, coords, sam)

        return aln

    def _parse_sam(self, sam, layer_names):
        aln = self.sam_to_aln(sam)

        layers = aln.to_pandas(layer_names, index=["seq.fwd", "seq.pos", "aln.id"])
        layers["track.id"] = self.in_id
        layers = layers.set_index("track.id", append=True)\
                       .reorder_levels(["track.id", "seq.fwd", "seq.pos", "aln.id"])

        attr = aln.attrs()
        aln_df = pd.DataFrame([attr], columns=attr._fields)
        aln_df["track.id"] = self.in_id
        aln_df.set_index(["track.id","id"], inplace=True)

        return aln_df, layers

    def query(self, layers, coords, index=["aln.id","seq.pos"], read_id=None, full_overlap=False):
        itr = self.input.fetch(coords.ref_name, coords.refs.min(), coords.refs.max())

        track_alns = defaultdict(list)
        track_layers = defaultdict(list)

        if read_id is not None and len(read_id) > 0:
            read_ids = set(read_id)
        else:
            read_ids = None

        for sam in itr:
            if read_ids is not None and sam.query_name is not read_ids:
                continue
            aln = self.sam_to_aln(sam)

            track_alns[self.in_id].append(aln.attrs())
            l = aln.to_pandas(layers, index=index)
            track_layers[self.in_id].append(l)

        track_alns = {
            t : pd.DataFrame(l, columns=Alignment.Attrs._fields)\
                  .set_index("id")
            for t,l in track_alns.items()}
        track_layers = {t : pd.concat(l) for t,l in track_layers.items()}
            
        return track_alns, track_layers

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
            next_alns,next_layers = self._parse_sam(sam, layers)

            next_aln = next_alns.iloc[0]
            next_ref = next_aln["ref_name"]
            next_start = next_layers.index.get_level_values("seq.pac").min()

            #we're losing alignments somewhere. getting out of sync
            #maybe right below
            #but also, not sure if I'm syncing inputs in "Tracks"
            #probably need to always advance the minimum ref pos/ID

            ret_layer_count = 0

            if prev_ref is not None and (next_ref != prev_ref or next_start - prev_start > 0) and len(new_alns) > self.prms.aln_chunksize:
                #pac = self.tracks.index.pos_to_pac(a["ref_name"], prev_start)

                aln_df = pd.concat([aln_df] + new_alns)#.sort_index()
                layer_df = pd.concat([layer_df] + new_layers).sort_index()
                new_alns = list()
                new_layers = list()
                
                ret_layers = layer_df.loc[slice(None),slice(None),prev_start:next_start-1]
                ret_alns = aln_df.loc[ret_layers.index.droplevel(["seq.fwd", "seq.pac"]).unique()]

                layer_df = layer_df.loc[slice(None),slice(None),next_start:]
                aln_df = aln_df.loc[layer_df.index.droplevel(["seq.fwd", "seq.pac"]).unique()]
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
        ret_alns = aln_df.loc[ret_layers.index.droplevel(["seq.fwd", "seq.pac"]).unique()]

        if len(ret_alns) > 0:
            yield (ret_alns, ret_layers)



    def close(self):
        if self.input is not None:
            self.input.close()
        if self.output is not None:
            self.output.close()

