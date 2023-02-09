
import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections
from collections import defaultdict
import scipy

from .io import SQL, TSV, Eventalign, Tombo, BAM, ModelTrainer, INPUTS, OUTPUTS, INPUT_PARAMS, OUTPUT_PARAMS
from .aln_track import AlnTrack
from .layers import LAYER_META, parse_layers
from ..index import load_index, RefCoord, str_to_coord
from ..pore_model import PoreModel
from ..read_index import ReadIndex, Fast5Reader
#from ..fast5 import Fast5Reader
from .. import config
from . import Bcaln

from time import time


_REFSTAT_AGGS = {
    "mean" : np.mean, 
    "median" : np.median, 
    "q5" : (lambda x: np.quantile(x, 0.05)),
    "q95" : (lambda x: np.quantile(x, 0.95)),
    "q25" : (lambda x: np.quantile(x, 0.25)),
    "q75" : (lambda x: np.quantile(x, 0.75)),
    "stdv" : np.std, 
    "var"  : np.var,
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : np.min, 
    "max"  : np.min,
}

LAYER_REFSTATS = {"min", "max", "mean", "median", "stdv", "var", "skew", "kurt", "q25", "q75", "q5", "q95"}
COMPARE_REFSTATS = {"ks"}
ALL_REFSTATS = LAYER_REFSTATS | COMPARE_REFSTATS

REFSTAT_LABELS = {
    "min" : "Minimum", 
    "max" : "Maximum", 
    "mean" : "Mean", 
    "median" : "Median", 
    "stdv" : "Std. Dev.", 
    "var" : "Variance", 
    "skew" : "Skewness", 
    "kurt" : "Kurtosis",
    "ks" : "KS",
    "q5" : "5% Quantile",
    "q25" : "25% Quantile",
    "q75" : "75% Quantile",
    "q95" : "95% Quantile",
}

BUILTIN_TRACKS = {"_refstats", "_layerstats", "_readstats"}

CMP_GROUPS = {"cmp", "bc_cmp"}

class RefstatsSplit:
    def __init__(self, stats, track_count):
        self.layer = [s for s in stats if s in LAYER_REFSTATS]
        self.compare = [s for s in stats if s in COMPARE_REFSTATS]

        self.layer_agg = [_REFSTAT_AGGS[s] for s in self.layer]

        if len(self.layer) + len(self.compare) != len(stats):
            bad_stats = [s for s in stats if s not in ALL_REFSTATS]
            raise ValueError("Unknown stats: %s (options: %s)" % (", ".join(bad_stats), ", ".join(ALL_REFSTATS)))

        if len(self.compare) > 0 and track_count != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(self.compare))

class Tracks:
    def __init__(self, *args, **kwargs):
        self.read_index = kwargs.get("read_index", None)
        if self.read_index is not None:
            del kwargs["read_index"]

        if len(args) > 0 and isinstance(args[0], Tracks):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

    def _init_new(self, *args, **kwargs):

        self.conf, self.prms = config._init_group("tracks", copy_conf=True, *args, **kwargs)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers, add_deps=False))

        self.alns = list()
        self._tracks = dict()
        self.new_alignment = False
        self.new_layers = set()

        
        if self.read_index is None:
            self.read_index = ReadIndex(read_filter=self.prms.read_filter)

        #def __init__(self, index_filename=None, file_paths=None, read_filter=None, file_suffix=".fast5"):

        self._fast5s = None

        self._init_io()

        self.set_layers(self.prms.layers)

        if self.prms.index_prefix is None:
            raise RuntimeError("Failed to load reference index")

        self._aln_track_ids = [t.id for t in self.alns]

        self.index = load_index(self.model.K, self.prms.index_prefix, kmer_shift=self.model.PRMS.shift)

        self.coords = self._ref_bounds_to_coords(self.prms.ref_bounds)

        #if self.coords is not None and  len(self._aln_track_ids) > 0:
        #    self.load()


    @property
    def fast5s(self):
        if self._fast5s is not None:
            return self._fast5s

        #for io in self.inputs:
        #    if isinstance(io, SQL):
        #        fast5_reads = list()
        #        fast5_reads.append(io.get_fast5_index(self._aln_track_ids))
        #        fast5_reads = pd.concat(fast5_reads)
        #        files = fast5_reads["filename"].unique()
        #        self._fast5s = Fast5Reader(
        #            index=fast5_reads, 
        #            conf=self.conf)

        #if self._fast5s is None:
        if len(self.conf.fast5_reader.fast5_files) > 0:
            self._fast5s = Fast5Reader(self.read_index, conf=self.conf)
        else: 
            return None

        for track in self.alns:
            track.fast5s = self._fast5s

        return self._fast5s

    def _init_slice(self, parent, coords, tracks):
        self.conf = parent.conf 
        self.prms = parent.prms

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers, add_deps=False))

        self.index = parent.index
        self.read_index = parent.read_index
        self._fast5s = parent._fast5s
        self._aln_track_ids = parent._aln_track_ids
        self.refstats = None
        self.new_alignment = parent.new_alignment
        self.new_layers = parent.new_layers

        self.inputs = parent.inputs
        self.output = parent.output
        self.output_track = parent.output_track

        self.coords = coords
        self._tracks = dict()
        self.alns = list()

        if tracks is not None:
            self._add_tracks(tracks)

    def _add_tracks(self, tracks):
        for name,track in tracks.items():
            self._add_track(name, track)

    def _add_track(self, name, track):
        if name in self._tracks:
            raise KeyError(f"Duplicate track name: {name}")
        self._tracks[name] = track
        
        if name in BUILTIN_TRACKS:
            setattr(self, name[1:], track)
        elif isinstance(track, AlnTrack):
            self.alns.append(track)
        else:
            raise ValueError("Unrecognized track type: " + str(track))


    @property
    def all_empty(self):
        for t in self.alns:
            if not t.empty:
                return False
        return True

    @property
    def any_empty(self):
        for t in self.alns:
            if t.empty:
                return True
        return False

    @property
    def track_names(self):
        return list(self._tracks.keys())

    def _ref_bounds_to_coords(self, ref_bounds):
        if ref_bounds is not None:
            if isinstance(ref_bounds, str):
                ref_bounds = str_to_coord(ref_bounds)
            elif isinstance(ref_bounds, tuple):
                ref_bounds = RefCoord(*ref_bounds)
            return self.index.get_coord_space(ref_bounds, self.conf.is_rna)
        return None

    def set_layers(self, layers):
        self.prms.layers = layers

        self.db_layers = list()
        self.fn_layers = list()
        self.cmp_layers = list()
        for group, layer in parse_layers(layers):
            if LAYER_META.loc[(group,layer),"fn"] is None:
                if group in CMP_GROUPS:
                    self.cmp_layers.append((group, layer))
                else:
                    self.db_layers.append((group, layer))
            else:
                self.fn_layers.append((group, layer))

    def __len__(self):
        return len(self._tracks)

    def keys(self):
        return self._tracks.keys()

    def __getitem__(self, i):
        return self._tracks[i]

    def _init_io(self):

        #tracks = list()

        self.inputs = list()
        self.bam_in = None

        track_count = 0

        for name,Cls in INPUTS.items():
            fnames = getattr(self.prms.io, name, None)
            if fnames is None: continue
            if isinstance(fnames, str):
                fnames = [fnames]
            assert isinstance(fnames, list)
                
            for filename in fnames:
                io = Cls(filename, False, self, track_count)

                p = io.conf.fast5_reader
                self.read_index.load_paths(p.fast5_files, p.recursive)
                self.read_index.load_index_file(p.fast5_index)

                self.conf.load_config(io.conf)
                self.inputs.append(io)

                if isinstance(io, BAM) and self.bam_in is None:
                    self.bam_in = io

                track_count += len(io.aln_tracks)
                for track in io.aln_tracks:
                    self._add_track(track.name, track)

        in_prms = [getattr(self.prms.io, p) is not None for p in INPUT_PARAMS]
        out_prms = [getattr(self.prms.io, p) is not None for p in OUTPUT_PARAMS]
        #if np.sum(in_prms) > 1:
        #    raise ValueError("No more than one input can be specified")
        if np.sum(out_prms) > 1:
            raise ValueError("No more than one output can be specified")

        if np.any(out_prms):
            out_format = OUTPUT_PARAMS[out_prms][0]
            filename = getattr(self.prms.io, out_format)
            new_track = True
            if out_format == "sql_out":
                self.output = SQL(filename, True, self, track_count)
            elif out_format == "tsv_out":
                self.output = TSV(filename, True, self, track_count)
            elif out_format == "eventalign_out":
                self.output = Eventalign(filename, True, self, track_count)
            elif out_format == "bam_out":
                self.output = BAM(filename, True, self, track_count)
            elif out_format == "model_dir":
                self.output = ModelTrainer(filename, True, self, track_count)

            track_count += len(self.output.aln_tracks)
            for track in self.output.aln_tracks:
                self._add_track(track.name, track)

                #tracks.append(self.output.aln_tracks)

            self.output_track = self.output.aln_tracks[0].name
            #for track in self.output.tracks:
            #    self.output_tracks[track.name] = track
        else:
            self.output = None
            self.output_track = self.inputs[0].aln_tracks[0].name

        self.model = None

        #for _,row in pd.concat(tracks).iterrows():
        #    conf = config.Config(toml=row["config"])
        #    self.conf.load_config(conf)
        #    track = AlnTrack(row["id"], row["name"], row["desc"], conf)
        #    self._add_track(track.name, track)

            #TODO each IO should have reference to its AlnTracks?
            #probably construct here in Tracks? or maybe in IO!
            #iterate through each input, it will populate its own AlnTracks

        for track in self.alns:
            if self.model is None:
                self.model = track.model
            elif self.model.K != track.model.K:
                raise ValueError("Cannot load tracks with multiple k-mer lengths (found K={self.model.K} and K={track.model.K}")

    def get_read_fast5(self, read_id):
        if self.read_index is not None and self.read_index.read_files is not None:
            return self.read_index.read_files.loc[read_id, "filename"]
        for io in self.inputs:
            if getattr(io, "read_id_in", None) == read_id:
                return io.fast5_in
        raise RuntimeError("Could not determine fast5 filename, may need to provide fast5 index (-x)")

    def aln_layers(self, layer_filter=None):
        ret = pd.Index([])
        for track in self.alns:
            layers = track.layers.columns
            if layer_filter is not None:
                layers = layers.intersection(layer_filter)
            ret = ret.union(layers)
        return ret
    
    def _track_or_default(self, track_name):
        if track_name is None:
            return self._tracks[self.output_track]
        elif track_name in self._tracks:
            return self._tracks[track_name]
        raise ValueError(f"Unknown track: {track_name}")

    def collapse_events(self, dtw, read=None):
        dtw["cuml_mean"] = dtw["length"] * dtw["current"]

        grp = dtw.groupby(level=0)

        lengths = grp["length"].sum()


        dtw = pd.DataFrame({
            "start"  : grp["start"].min().astype("uint32"),
            "length" : lengths.astype("uint32"),
            "current"   : grp["cuml_mean"].sum() / lengths,
            "kmer" : grp["kmer"].first(),
            "events" : grp["start"].count(),
        })

        if read is not None:
            #dtw["stdv"] = np.std(read.signal[dtw["start"]:dtw["start"]+dtw["length"]], axis=1)
            dtw["stdv"] = [
                np.std(
                    read.get_norm_signal(dtw.loc[i, "start"], dtw.loc[i, "start"]+dtw.loc[i, "length"])
                ) for i in dtw.index
            ]

        skip_counts = dtw["start"].value_counts().loc[dtw["start"]].to_numpy()
        dtw["events"] /= skip_counts

        return dtw.sort_index()

    def write_dtw_events(self, events=None, track_name=None, aln_id=None, read=None):
        if np.any(events.index.duplicated()):
            events = self.collapse_events(events, read=read)

        self.add_layers("dtw", events, track_name, aln_id, False, read)
        
    def add_layers(self, group, layers, track_name=None, aln_id=None, overwrite=False, read=None):
        track = self._track_or_default(track_name)

        if layers.index.names[0] == "mref":
            layers = layers.set_index(self.index.mref_to_ref(layers.index))
        elif layers.index.names[0] == "pac":
            layers = layers.set_index(self.index.pac_to_ref(layers.index))

        layers.index.names = ("ref",)

        track.add_layer_group(group, layers, aln_id, overwrite)

        self.new_layers.add(group)


    def write_alignment(self, track_name=None):
        track = self._track_or_default(track_name)

        if self.output is not None:
            out = self.output
        else:
            out = self.inputs[0]
            #raise RuntimeError("Error writing output")

        if self.new_alignment:
            out.write_alignment(track.alignments)

        if len(self.new_layers) > 0 and len(track.layers) > 0:
            out.write_layers(track, self.new_layers)

    def set_read(self, read):
        self.output.read = read

    def init_alignment(self, read_id, fast5, coords, layers={}, read=None, bam=None, track_name=None):
        track = self._track_or_default(track_name)
        self.new_alignment = True
        self.new_layers = set()

        aln_id = self.output.init_alignment(read_id, fast5, read=read, bam=bam)

        if len(layers) > 0:
            starts = None
            for group,vals in layers.items():
                if "start" in vals.columns:
                    starts = vals["start"]
                    break

            if starts is None:
                raise ValueError("Must initialize alignment from DataFrame with start column")

            samp_start = starts.min()
            samp_end = starts.max()
        else:
            samp_start = samp_end = None

        track.alignments = pd.DataFrame({
                "id" : [aln_id],
                "track_id" : [track.id],
                "read_id" : [read_id],
                "ref_name" : [coords.ref_name],
                "ref_start" : [coords.refs.start],
                "ref_end" : [coords.refs.stop],
                "fwd" :     [coords.fwd],
                "samp_start" : [samp_start],
                "samp_end" : [samp_end],
                "tags" : [""]}).set_index("id")
        

        track.layers = None

        track.coords = coords

        for group,vals in layers.items():
            if group == "dtw":
                self.write_dtw_events(vals, track.name, aln_id)
            else:
                self.add_layers(group, vals, track.name)

        return aln_id, coords

    @property
    def input_count(self):
        return len(self.alns)

    def _verify_read(self):
        if len(self._aln_track_ids) == 0:
            raise RuntimeError("No input tracks have been loaded")

    #TODO read_ids, track_alns, max_cov(?)
    def slice(self, ref_start=None, ref_end=None, reads=None, order=["fwd","ref_start"], tracks=None, full_overlap=False, shared_reads=False):
        if self.coords is None:
            raise IndexError("Cannot slice empty Tracks")

        hasbounds = (ref_start is not None, ref_end is not None)
        if np.all(hasbounds):
            coords = self.coords.ref_slice(ref_start, ref_end)
        elif np.any(hasbounds):
            raise IndexError(f"Invalid bounds {ref_start}-{ref_end}")
        else:
            coords = self.coords

        stranded = True
        fwd = None

        if full_overlap:
            reads = self.get_full_overlap(reads)

        if shared_reads:
            if reads is None:
                reads = self.get_shared_reads()
            else:
                reads = self.get_shared_reads().intersection(reads)
            order = "read_id"
                
        if tracks is None:
            track_names = set(self.track_names)
        else:
            track_names = set(tracks)
        tracks = dict()
        for name,track in self._tracks.items():
            if name not in track_names: continue
            if isinstance(track, pd.DataFrame):
                tracks[name] = track.loc[coords.refs]

            elif isinstance(track, AlnTrack):
                tracks[name] = track.slice(coords, reads=reads, order=order)

        return Tracks(self, coords, tracks)

    def get_shared_reads(self):
        read_ids = pd.Index(self.alns[0].read_ids)
        for track in self.alns[1:]:
            read_ids = read_ids.intersection(track.read_ids)
        return read_ids

    def get_all_reads(self):
        read_ids = pd.Index(self.alns[0].read_ids)
        for track in self.alns[1:]:
            if not track.empty:
                read_ids = read_ids.union(track.read_ids)
        return read_ids

    def get_full_overlap(self, read_ids=None):
        rmin = self.coords.refs.min()
        rmax = self.coords.refs.max()
        for track in self.alns:
            alns = track.alignments[(track.alignments["ref_start"] <= rmin) & (track.alignments["ref_end"] >= rmax)]
            if read_ids is None:
                read_ids = pd.Index(alns["read_id"])
            else:
                read_ids = read_ids.intersection(alns["read_id"])
        return read_ids
            
    def load(self, ref_bounds=None, full_overlap=None, read_filter=None, load_mat=False):
        self._verify_read()
        self.new_alignment = True
        self.new_layers = set()

        if ref_bounds is not None:
            self.coords = self._ref_bounds_to_coords(ref_bounds)

        if read_filter is None:
            read_filter = self.prms.read_filter

        if self.coords is None:
            raise ValueError("Must set ref bounds")

        if full_overlap is None:
            full_overlap = self.prms.full_overlap

        if load_mat is None:
            load_mat = self.prms.load_mat

        for io in self.inputs:
            alignments, layers = io.query(self.db_layers, self._aln_track_ids, self.coords, full_overlap=full_overlap, read_id=read_filter)

            for track in io.aln_tracks:
                track_alns = alignments.loc[track.id]
                track_layers = layers.loc[track.id].droplevel(0)

                track.set_data(self.coords, track_alns, track_layers)
                track.calc_layers(self.fn_layers)
            
        #TODO get this back
        self.load_compare(alignments.index.droplevel(0).to_numpy())
        self.calc_refstats()

        for track in self.alns:
            track.calc_layers(self.fn_layers)
            if load_mat:
                track.load_mat()

        return self.alns

    def load_compare(self, aln_ids=None):
        if len(self.cmp_layers) == 0:
            return

        #TODO handle multiple inputs properly
        io = self.inputs[0]

        self.cmp = io.query_compare(self.cmp_layers, self._aln_track_ids, self.coords, aln_ids)


        if self.cmp is None:
            return

        groups = self.cmp.index.get_level_values("group_b").unique()
        if "bcaln" in groups:
            bcalns = self.cmp.loc[(slice(None), slice(None), slice(None), "bcaln"),:]
        else:
            bcalns = None
        if "dtw" in groups:
            dtws = self.cmp.loc[(slice(None), slice(None), slice(None), "dtw"),:]
        else:
            dtws = None

        for track in self.alns:
            def _add_group(group, df):
                df = df.reset_index(["aln_b", "group_b"])
                df = df[df.index.get_level_values("pac").isin(track.layer_pacs)]
                df.rename(index=track.coords.pac_to_ref, level=0, inplace=True)
                df.index.names = ["ref", "aln_id"]
                df = pd.concat({group : df.reindex(track.layers.index)}, axis=1)
                track.layers = pd.concat([track.layers, df], axis=1).dropna(axis=1,how="all")

            #try:
            if bcalns is not None:
                _add_group("bc_cmp", bcalns)
            if dtws is not None:
                _add_group("cmp", dtws)
            #except:
            #    sys.stderr.write("Failed to write compare group\n")
            #    sys.stderr.write(str(track.alignments))

    def calc_bcaln(self, bam, read, track_name=None):
        track = self._track_or_default(track_name)

        read_id = bam.query_name

        #TODO move to tracks.init_bcaln (maybe rename bcaln)
        #construct from BAM (with move table, or FAST5)
        #need to standardize ref_gaps, maybe event_gaps too
        #ideally in generalized C++ CoordBounds or whatever

        bcaln = Bcaln(self.conf, self.index, read, bam, self.coords)

        if bcaln.empty:
            return None, None

        if track.empty or not read_id in track.alignments["read_id"]:
            aln_id, aln_coords = self.init_alignment(read_id, read.filename, bcaln.coords, {"bcaln" : bcaln.df}, bam=bam) #, read=signal
        else:
            self.add_layers("bcaln", bcaln.df, track.name, aln_id)
        
        return bcaln, aln_coords

    def calc_compare(self, group_b, single_track, save):
        if len(self.alns) > 0:
            alns = self.alns
        #elif len(self.output_tracks) > 0:
        #    alns = list(self.output_tracks.values())
        else:
            raise ValueError("Must input at least one track")

        #if single_track:
        #    if self.output_track is None:
        #        track_a = self.alns[0]
        #    else:
        #        track_a = self._tracks[self.output_track]
        #    track_b = track_a
        #elif len(self.alns) == 1:
        #    track_a = self.alns[0]
        #elif self.output_track is not None:
        
        if self.output_track is not None:
            track_a = self._tracks[self.output_track]
        elif single_track or len(self.alns) == 1:
            track_a = self.alns[0]
        else:
            track_a = self.alns[0]

        if single_track or len(self.alns) == 1:
            track_b = track_a
        else:
            for track_b in self.alns:
                if track_b != track_a: break

        cols = track_a.layers.columns.get_level_values("group").unique()
        if (group_b == "dtw" and "cmp" in cols) or (group_b == "bcaln" and "bc_cmp" in cols):
            sys.stderr.write(f"Read already has compare group. Skipping\n")
            return None

        if group_b == "dtw":
            if track_a == track_b:
                raise ValueError("Must input exactly two tracks to compare dtw alignments")

            df = track_a.cmp(track_b, True, True)

        elif group_b == "bcaln":
            df = track_a.bc_cmp(track_b, True, True)

        df = df.dropna(how="all")
        self.add_layers("cmp", df)

    def calc_refstats(self, cov=False):
        if self.prms.refstats is None or len(self.prms.refstats) == 0 or len(self.prms.refstats_layers) == 0 or self.all_empty:
            self.refstats = None
            return None

        stats = RefstatsSplit(self.prms.refstats, len(self.alns))

        refstats = dict()
        grouped = [
            t.layers[self.prms.refstats_layers].groupby(level=0)
            for t in self.alns]

        for track,groups in zip(self.alns, grouped):
            if track.empty:
                refstats[track.name] = None
                continue

            refstats[track.name] = groups.agg(stats.layer_agg)
            rename = ({
                old[-1] : new
                for old,new in zip(refstats[track.name].columns, stats.layer)
            })
            refstats[track.name].rename(columns=rename, inplace=True)

            if cov:
                refstats[track.name].insert(0, "cov", groups.size())

        if len(stats.compare) > 0:
            if self.any_empty:
                refstats["ks"] = None
            
            else:
                groups_a, groups_b = grouped
                refs_a = self.alns[0].layers.index.unique("ref")
                refs_b = self.alns[1].layers.index.unique("ref")
                refs = refs_a.intersection(refs_b)
                cmps = {l : defaultdict(list) for l in self.prms.refstats_layers}
                for ref in refs:
                    track_a = groups_a.get_group(ref)
                    track_b = groups_b.get_group(ref)
                    for layer in self.prms.refstats_layers:
                        a = track_a[layer]
                        b = track_b[layer]
                        for stat in stats.compare:
                            ks = scipy.stats.stats.ks_2samp(a,b,mode="asymp")
                            cmps[layer]["stat"].append(ks.statistic)
                            cmps[layer]["pval"].append(ks.pvalue)

                refstats["ks"] = pd.concat({k : pd.DataFrame(index=refs, data=c) for k,c in cmps.items()}, axis=1) 
        
        if np.any([df is None for df in refstats.values()]):
            columns = None
            for df in refstats.values():
                if df is not None:
                    columns = df.columns
                    break
            for name,df in refstats.items():
                if df is None:
                    refstats[name] = pd.DataFrame(columns=columns)
        refstats = pd.concat(refstats, axis=1, names=["track", "group", "layer", "stat"])

        #TODO make verbose ref indexing
        #refstats.index = self.alns[0].coords.mref_to_ref_index(refstats.index, multi=verbose_refs)

        self.refstats = refstats.dropna()

        self._tracks["_refstats"] = self.refstats

        return refstats

    def iter_refs(self, ref_bounds=None):
        if ref_bounds is not None:
            coords = self._ref_bounds_to_coords(ref_bounds)
        else:
            coords = self.coords

        def get_full_coords(pac):
            rid = self.index.pac_to_ref_id(pac)
            ref_len = self.index.get_ref_len(rid)
            ref_name = self.index.get_ref_name(rid)
            seq_refs = RefCoord(ref_name, 0, ref_len)
            return self.index.get_coord_space(seq_refs, self.conf.is_rna, load_kmers=False)

        def next_coords(seq_coords, pac, fwd, chunk_hasnext):
            if len(pac) > 1 and chunk_hasnext:
                pac = pac[:-1]

            idx = pd.MultiIndex.from_product([[int(fwd)], pac])

            coords = seq_coords.pac_intersect(idx)
            if coords is None:
                seq_coords = get_full_coords(pac[0])
                coords = seq_coords.pac_intersect(idx)

            return seq_coords, coords

        #NEW CODE START 

        layer_iters = [
            io.iter_refs(
                self.db_layers, 
                self._aln_track_ids, 
                coords=coords, 
                chunksize=self.prms.io.ref_chunksize)
            for io in self.inputs]

        chunks = [next(i) for i in layer_iters]
        chunk_hasnext = np.ones(len(chunks), bool)

        pac_min = min([l.index.get_level_values("pac")[0] for a,l in chunks])
        seq_coords = get_full_coords(pac_min)


        while np.any([len(chunk[0]) > 0 for chunk in chunks]):
            pac_start = np.inf
            pac_end = np.inf

            all_pacs = pd.Index([])

            strand = -1
            all_fwd = True
            all_rev = True

            #TODO simplify! assume inputs segment by ref coord
            for i in range(len(chunks)):
                chunk_alignments, chunk_layers = chunks[i]
                last_chunk = False

                pacs = chunk_layers.index.get_level_values("pac").unique()
                while len(pacs) <= 1:
                    alns, layers = next(layer_iters[i], (None, None))
                    if alns is None:
                        chunk_hasnext[i] = False
                        break

                    chunk_alignments = pd.concat([chunks[i][0], alns])
                    chunk_alignments = chunk_alignments[~chunk_alignments.index.duplicated()]
                    chunk_layers = pd.concat([chunk_layers, layers])
                    chunks[i] = (chunk_alignments,chunk_layers)
                    pacs = chunk_layers.index.get_level_values("pac").unique()

                idx = chunk_layers.index
                if all_fwd or all_rev:
                    fwds = idx.get_level_values("fwd").unique()
                    all_fwd = all_fwd and np.all(fwds == 1)
                    all_rev = all_rev and np.all(fwds == 0)
                pacs = idx.get_level_values("pac").unique()

                all_pacs = all_pacs.union(pacs)
                pac_start = min(pac_start, pacs.min())

                #TODO this part is still key. only output up until the end of the least full chunk
                pac_end = min(pac_end, pacs.max())+1

            chunk_pacs = pd.RangeIndex(pac_start, pac_end)
            
            #chunk_pacs = all_pacs[all_pacs <= pac_end].sort_values()

            if not (all_fwd or all_rev):
                strands = [False, True]
                fwd = None
            else:
                strands = [all_fwd]
                fwd = all_fwd


            for fwd in strands:

                seq_coords, coords = next_coords(seq_coords, chunk_pacs, fwd, chunk_hasnext[i])
                coords.set_kmers(self.index.mrefs_to_kmers(coords.mrefs, self.conf.is_rna, False))

                #TODO probably don't need masks anymore
                masks = [
                    l.index.get_level_values("pac").isin(coords.pacs) & (l.index.get_level_values("fwd") == fwd)
                    for a,l in chunks]


                layers = pd.concat([l[mask] for (a,l),mask in zip(chunks, masks)])

                alns = pd.concat([a.loc[l[mask].index.droplevel(["fwd","pac"]).unique()] for (a,l),mask in zip(chunks, masks)])
                aln_ids = alns.index

                ret = self._tables_to_tracks(coords, alns, layers)

                #leftover collection
                for i in range(len(chunks)):
                    chunk_alns, chunk_layers = chunks[i]

                    #if last_chunk:
                    #    chunk_alns = chunk_alns.iloc[:0]
                    #    chunk_layers = chunk_layers.iloc[:0]
                    #else:
                    chunk_layers = chunk_layers[~masks[i]]
                    if len(chunk_layers) > 0:
                        ids = chunk_layers.index.droplevel(["fwd","pac"]).unique()
                        chunk_alns = chunk_alns.loc[ids]
                    else:
                        chunk_alns = chunk_alns.iloc[:0]

                    chunks[i] = (chunk_alns,chunk_layers)

                yield ret

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None):

        if read_filter is None and self.read_index is not None:
            read_filter = self.read_index.read_filter
        
        if ref_bounds is not None and not isinstance(ref_bounds, RefCoord):
            ref_bounds = RefCoord(ref_bounds)

        if (self.coords is None or
            (read_filter is not None and 
             len(self.get_all_reads().intersection(read_filter)) < len(read_filter)) or
            (ref_bounds is not None and not self.coords.contains(ref_bounds))):
            gen = self.iter_reads_db(read_filter, ref_bounds, full_overlap, max_reads)
        else:
            gen = self.iter_reads_slice(read_filter, ref_bounds)

        for read_id,chunk in gen:
            #print(read_id, [len(a.alignments) for a in chunk.alns])
            yield read_id,chunk
            
    def iter_reads_slice(self, reads=None, ref_bounds=None):
        all_reads = self.get_all_reads()
        if reads is not None:
            all_reads = all_reads.intersection(reads)

        if ref_bounds is None:
            ref_start = ref_end = None
        else:
            ref_start = ref_bounds.start
            ref_end = ref_bounds.end
        
        for read_id in all_reads:
            yield read_id, self.slice(ref_start, ref_end, [read_id])

    def iter_reads_db(self, reads, ref_bounds, full_overlap, max_reads):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if reads is None:
            reads = self.prms.read_filter
        if max_reads is None:
            max_reads = self.prms.max_reads

        #TODO handle multiple inputs properly
        for io in self.inputs:
            #aln_iter = self.input.iter_alns(
            aln_iter = io.iter_alns(
                self.db_layers, 
                self._aln_track_ids,
                self.coords,
                read_id=reads,
                full_overlap=full_overlap,
                ref_index=self.index)

            for alignments,layers in aln_iter:
                #print("aln", alignments)
                #print("layer", layers.index.unique("aln_id"))
                for ref_name,ref_alns in alignments.groupby("ref_name"):
                    coords = self._alns_to_coords(ref_alns)
                    cache = self._tables_to_tracks(coords, ref_alns, layers)
                    for ret in cache.iter_reads_slice():
                        yield ret
        
    def _tables_to_tracks(self, coords, alignments, layers):
        tracks = dict()
        self.new_alignment = False
        self.new_layers = set()

        layer_alns = layers.index.get_level_values("aln_id")

        if self.prms.shared_refs_only or self.prms.min_coverage > 1:
            #track_covs = alignments.loc[layer_alns, ["track_id"]] 
            #track_covs = track_covs.set_index(layers.index)
            #track_covs = track_covs.reset_index("aln_id", drop=True) 
            #track_covs = track_covs.set_index("track_id", append=True) 
            track_covs = layers.index.droplevel("aln_id")
            track_covs = track_covs.value_counts()

            mask = track_covs >= self.prms.min_coverage
            if not np.any(mask):
                idx = layers.index[:0]
            elif self.prms.shared_refs_only:
                track_counts = pd.MultiIndex.from_tuples(track_covs[mask].index) \
                                   .droplevel(0) \
                                   .value_counts()

                shared = track_counts.index[track_counts == len(self.alns)]
                if len(shared) == 0: 
                    idx = layers.index[:0]
                else:
                    idx = pd.MultiIndex.from_tuples(shared)
            else:
                idx = track_covs[mask].index.unique()

            layers = layers.loc[(slice(None),idx.get_level_values(0),idx.get_level_values(1),slice(None))]
            layer_alns = layers.index.droplevel(["fwd","pac"])
            alignments = alignments.loc[layer_alns.unique()]

        aln_groups = alignments.index.unique("track_id")
        for parent in self.alns:
            if parent.id in aln_groups:
                track_alns = alignments.loc[parent.id]
                track_layers = layers.loc[parent.id].droplevel("fwd")
            else:
                track_alns = alignments.iloc[:0] 
                track_layers = layers.iloc[:0]   

            track = AlnTrack(parent, coords, track_alns, track_layers)

            #if not track.empty:
            track.calc_layers(self.fn_layers)

            tracks[parent.name] = track

        tracks = Tracks(self, coords, tracks)

        if not tracks.all_empty:
            tracks.load_compare(alignments.index.to_numpy())
        tracks.calc_refstats()

        return tracks

    def _alns_to_coords(self, alns):
        ref_coord = RefCoord(
            alns["ref_name"].iloc[0],
            alns["ref_start"].min(),
            alns["ref_end"].max())
        return self.index.get_coord_space(
            ref_coord, self.conf.is_rna, load_kmers=True, kmer_trim=True)

    
    def close(self):
        for io in self.inputs:
            io.close()
        if self.output is not None and (len(self.inputs) == 0 or self.output != self.inputs[0]):
            self.output.close()
