
import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections
from collections import defaultdict, namedtuple
import scipy

from .io import SQL, TSV, Eventalign, Tombo, BAM, ModelTrainer, INPUTS, OUTPUTS, INPUT_PARAMS, OUTPUT_PARAMS
from .aln_track import AlnTrack
#from .layers import LAYER_META#, parse_layers
from ..aln import Alignment, Sequence, AlnDF, parse_layers, LAYER_META
from ..index import load_index, RefCoord, str_to_coord
from ..pore_model import PoreModel, WORKFLOW_PRESETS
from ..read_index import ReadIndex, Fast5Reader, Slow5Reader
#from ..fast5 import Fast5Reader
from .. import config

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

_REFSTAT_CMPS = {
    "ks" : lambda df: scipy.stats.stats.ks_2samp(df.loc[df.index.unique(0)],df.loc[df.index.unique(1)],mode="asymp")
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

CMP_GROUPS = {"cmp", "mvcmp"}

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

    def _init_new(self, *args, model=None, **kwargs):

        self.conf, self.prms = config._init_group("tracks", copy_conf=True, *args, **kwargs)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers))

        self.alignments = None
        self.layers = None

        self.new_alignment = False
        self.new_layers = set()

        self.alns = list()
        self.models = dict()
        self._tracks = dict()
        self.new_alignment = False
        self.new_layers = set()
        
        if self.read_index is None:
            self.read_index = ReadIndex(self.conf.read_index.paths, self.prms.read_filter, self.conf.read_index.read_index, self.conf.read_index.recursive)

        if model is not None:
            self.set_model(model)
        else:
            pm = self.conf.pore_model
            if not pm.has_workflow():
                flowcell, kit = self.read_index.default_model
                if flowcell is not None and kit is not None:
                    pm.flowcell = flowcell
                    pm.kit = kit

            if len(pm.name) == 0 and pm.has_workflow():
                pm.name = PoreModel.PRESET_MAP.loc[pm.get_workflow(), "preset_model"]
            self.model = None

        self._init_io()

        self.set_layers(self.prms.layers)

        if self.prms.index_prefix is None:
            raise RuntimeError("Failed to load reference index")

        self._aln_track_ids = [t.id for t in self.alns]

        self.index = load_index(self.model, self.prms.index_prefix)
        self.refstats = None

        #self.coords = self._ref_bounds_to_coords(self.prms.ref_bounds)
        self.coords = self.prms.ref_bounds

        if self.coords is not None and  len(self._aln_track_ids) > 0:
            self.load()

    def _init_slice(self, parent, coords, reads=None):
        self.conf = parent.conf 
        self.prms = parent.prms

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers))

        self.index = parent.index
        self.read_index = parent.read_index
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
        self.model = parent.model

        self._add_tracks(parent._tracks)

        if parent.alignments is not None:
            mask = np.minimum(parent.alignments["ref_start"], coords.start) < np.maximum(parent.alignments["ref_end"], coords.end)
            if reads is not None:
                mask &= parent.alignments["read_id"].isin(reads)
            self.alignments = parent.alignments[mask] #pd.concat(track_alns, axis=0, names=["track", "id"])
            self.layers = parent.layers.reset_index(level="seq.pos").loc[self.alignments.index].set_index("seq.pos",append=True)
        else:
            self.alignments = None
            self.layers = None

        #.reorder_levels(["track.name","aln.id","seq.pos"]).loc[self.alignments.index].reorder_levels(["track.name","seq.pos","aln.id"]) #pd.concat(track_layers, axis=0, names=["track.name", "seq.pos", "aln.id"])
    
        if self.alignments is not None:
            self._init_mat()

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
        #elif isinstance(track, tuple):
            self.alns.append(track)
        else:
            raise ValueError("Unrecognized track type: " + str(track))

    @property
    def track_count(self):
        return len(self.alns)

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
            if group in CMP_GROUPS:
                self.cmp_layers.append((group, layer))
            else:
                self.db_layers.append((group, layer))

    def __len__(self):
        return len(self._tracks)

    def keys(self):
        return self._tracks.keys()

    def __getitem__(self, i):
        return self._tracks[i]

    def _init_io(self):

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

                p = io.conf.read_index
                if p.paths is not None:
                    self.read_index.load_paths(p.paths, p.recursive)
                self.read_index.load_index_file(p.read_index)

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

            self.output_track = self.output.aln_tracks[0].name
            #for track in self.output.tracks:
            #    self.output_tracks[track.name] = track
        else:
            self.output = None
            self.output_track = self.inputs[0].aln_tracks[0].name

            #TODO each IO should have reference to its AlnTracks?
            #probably construct here in Tracks? or maybe in IO!
            #iterate through each input, it will populate its own AlnTracks

        for track in self.alns:
            if self.model is None:
                self.model = track.model
            elif track.model is not None and self.model.K != track.model.K:
                raise ValueError("Cannot load tracks with multiple k-mer lengths (found K={self.model.K} and K={track.model.K}")

        #pm = self.conf.pore_model
        #has_flowkit = not (len(pm.flowcell) == 0 or len(pm.kit) == 0)
        #if not has_flowkit:
        #    flowcell, kit = self.read_index.default_model
        #    if not (kit is None or flowcell is None):
        #        if len(pm.flowcell) == 0:
        #            pm.flowcell = flowcell
        #        if len(pm.kit) == 0:
        #            pm.kit = kit
        #        has_flowkit = True
        ##PoreModel.PRESET_MAP.loc[(pm.flowcell, pm.kit), "preset_model"]

        #if self.model is None and has_flowkit:
        #    if len(self.conf.pore_model.name) == 0:
        #        self.model = PoreModel(PoreModel.PRESET_MAP.loc[(pm.flowcell, pm.kit), "preset_model"])
        #    else:
        #        self.model = PoreModel(params=self.conf.pore_model)
        #    for track in self.alns:
        #        track.model = self.model

    def set_model(self, model):
        self.conf.pore_model = model.PRMS
        self.model = model
        for track in self.alns:
            track.model = model

    def get_read_fast5(self, read_id):
        if self.read_index is not None and self.read_index.read_files is not None:
            return self.read_index.read_files.loc[read_id]
        for io in self.inputs:
            if getattr(io, "read_id_in", None) == read_id:
                return io.fast5_in
        raise RuntimeError("Could not determine fast5 filename, may need to provide fast5 index (-x)")

    def aln_layers(self, layer_filter=None):
        ret = self.layers.columns
        if layer_filter is None:
            return ret
        return ret.intersection(layer_filter)
        #ret = pd.Index([])
        #for track in self.alns:
        #    layers = track.layers.columns
        #    if layer_filter is not None:
        #        layers = layers.intersection(layer_filter)
        #    ret = ret.union(layers)
        #return ret
    
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
            dtw["current_sd"] = [
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

        if layers.index.names[0] == "mpos":
            layers = layers.set_index(self.index.mpos_to_pos(layers.index))
        elif layers.index.names[0] == "seq.pac":
            layers = layers.set_index(self.index.pac_to_pos(layers.index))

        layers.index.names = ("seq.pos",)

        track.add_layer_group(group, layers, aln_id, overwrite)

        self.new_layers.add(group)


    def init_alignment(self, track_name, aln_id, read, ref_id, coords, sam=None):
        #TODO make model name paramter, initialize with track specific model
        track = self._track_or_default(track_name)
        model = track.model.instance
        #seq = self.index.instance.get_kmers(model, ref_id, coords, self.conf.is_rna)
        seq = self.index.query(coords)
        seq = Sequence(seq, self.index.get_pac_offset(ref_id))
        return Alignment(aln_id, read, seq, sam)

    def write_alignment(self, aln):
        if self.output is not None:
            out = self.output
        else:
            out = self.inputs[0]

        out.write_alignment(aln)

    def set_read(self, read):
        self.output.read = read

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
            coords = self.coords.intersection(RefCoord(self.coords.name, ref_start, ref_end))
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
                tracks[name] = track #.slice(coords, reads=reads, order=order)

        return Tracks(self, coords, reads)

    def get_shared_reads(self):
        all_ids = self.alignments["read_id"]
        read_ids = pd.Index(all_ids.loc[self.alns[0].id])
        for track in self.alns[1:]:
            read_ids = read_ids.intersection(all_ids.loc[track.name])
        return read_ids

    def get_all_reads(self):
        all_ids = self.alignments["read_id"]
        read_ids = pd.Index(all_ids.loc[self.alns[0].name])
        for track in self.alns[1:]:
            read_ids = read_ids.union(all_ids.loc[track.name])
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

        track_alns = dict()
        track_layers = dict()

        for io in self.inputs:
            alns, layers = io.query(self.conf.tracks.layers, self.coords, ["aln.id","seq.pos"], full_overlap=full_overlap, read_id=read_filter)

            #for track in io.aln_tracks:
            #    track.set_data(self.coords, alns[track.name], layers[track.name])

            track_alns.update(alns)
            track_layers.update(layers)

        #TODO eliminate need for AlnTrack
        #just use self.layers, self.aln_layers
        #will need to initialize Alignment from these dataframes for Dotplot
        self.layers = pd.concat(track_layers, axis=0, names=["track.name", "aln.id", "seq.pos"])
        self.alignments = pd.concat(track_alns, axis=0, names=["track", "id"]).sort_index()#.sort_values(["fwd","ref_start"])
    
        self._init_mat()
            
        #self.load_compare(alignments.index.droplevel(0).to_numpy())
        self.calc_refstats()

        for track in self.alns:
            track.calc_layers(self.fn_layers)
            if load_mat:
                track.load_mat()

        return self.alns

    def _init_mat(self):
        df = self.layers.reset_index()

        self.mat = df.pivot(index=["track.name","aln.id"], columns=["seq.pos"]) 
        self.mat = self.mat.rename_axis(("group","layer","seq.pos"), axis=1) 
        self.mat = self.mat.reindex(pd.RangeIndex(self.coords.start, self.coords.end), axis=1, level=2) 
        self.mat = self.mat.sort_index(axis=1)

        order = self.alignments.sort_values(["fwd", "ref_start"]).index

        self.mat = self.mat.reindex(order, copy=False)

        self.width = len(self.coords)
        self.height = len(self.alignments)

    def load_compare(self, aln_ids=None):
        if len(self.cmp_layers) == 0:
            return

        #TODO handle multiple inputs properly
        io = self.inputs[0]

        self.cmp = io.query_compare(self.cmp_layers, self._aln_track_ids, self.coords, aln_ids)


        if self.cmp is None:
            return

        groups = self.cmp.index.get_level_values("group_b").unique()
        if "moves" in groups:
            movess = self.cmp.loc[(slice(None), slice(None), slice(None), "moves"),:]
        else:
            movess = None
        if "dtw" in groups:
            dtws = self.cmp.loc[(slice(None), slice(None), slice(None), "dtw"),:]
        else:
            dtws = None

        for track in self.alns:
            def _add_group(group, df):
                df = df.reset_index(["aln_b", "group_b"])
                df = df[df.index.get_level_values("seq.pac").isin(track.layer_pacs)]
                df.rename(index=track.coords.pac_to_pos, level=0, inplace=True)
                df.index.names = ["seq.pos", "aln.id"]
                df = pd.concat({group : df.reindex(track.layers.index)}, axis=1)
                track.layers = pd.concat([track.layers, df], axis=1).dropna(axis=1,how="all")

            #try:
            if movess is not None:
                _add_group("mvcmp", movess)
            if dtws is not None:
                _add_group("cmp", dtws)
            #except:
            #    sys.stderr.write("Failed to write compare group\n")
            #    sys.stderr.write(str(track.alignments))

    def calc_compare(self, group_b, single_track, save):
        if len(self.alns) > 0:
            alns = self.alns
        else:
            raise ValueError("Must input at least one track")
        
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
        if (group_b == "dtw" and "cmp" in cols) or (group_b == "moves" and "mvcmp" in cols):
            sys.stderr.write(f"Read already has compare group. Skipping\n")
            return None

        if group_b == "dtw":
            if track_a == track_b:
                raise ValueError("Must input exactly two tracks to compare dtw alignments")

            df = track_a.cmp(track_b, True, True)

        elif group_b == "moves":
            df = track_a.mvcmp(track_b, True, True)

        df = df.dropna(how="all")
        self.add_layers("cmp", df)

    def calc_refstats(self, cov=False):
        if self.prms.refstats is None or len(self.prms.refstats) == 0 or len(self.prms.refstats_layers) == 0 or self.alignments is None:
            self.refstats = None
            return None

        stats = RefstatsSplit(self.prms.refstats, len(self.alns))

        refstats = dict()
        #grouped = [
        #    t.layers[self.prms.refstats_layers].groupby(level=0)
        #    for t in self.alns]

        groups = self.layers[self.prms.refstats_layers].groupby(level=["track.name", "seq.pac", "seq.fwd"])

        refstats = groups.agg(stats.layer_agg).reset_index().pivot(index=["seq.pac","seq.fwd"], columns="track.name")
        #rename = ({
        #    old[-1] : new
        #    for old,new in zip(refstats.columns, stats.layer)
        #})
        #refstats.rename(columns=rename, inplace=True)
        if cov:
            refstats[track.name].insert(0, "cov", groups.size())

        if len(stats.compare) > 0:
            if len(self.layers.index.unique(0)) != 2:
                raise ValueError("Exactly two tracks must be input for comparison")
                #refstats["ks"] = None
            
            else:
                cmp_groups = self.layers[self.prms.refstats_layers].groupby(level=["seq.pac","seq.fwd"])
                def ks(df):
                    idx = df.index.levels[0]
                    d = dict()
                    for col in df.columns:
                        vals = scipy.stats.stats.ks_2samp(df.loc[idx[0],col],df.loc[idx[1],col],mode="asymp")
                        d[col+("ks",)] = pd.Series(vals, index=["stat","pval"], name=df.name)
                        #return pd.Series(v, index=pd.MultiIndex.from_product([["ks"],["stat","pval"]]), name=df.name)
                    return pd.concat(d, axis=0)

                df = cmp_groups.apply(ks)
                refstats = pd.concat([refstats, df], axis=1)

        self.refstats = refstats.dropna()

        self._tracks["_refstats"] = self.refstats

        return refstats

    def iter_refs(self, ref_bounds=None):
        layer_iters = [
            io.iter_refs(
                self.db_layers, 
                self._aln_track_ids, 
                #coords=coords, 
                chunksize=self.prms.io.ref_chunksize)
            for io in self.inputs]

        chunks = [next(i) for i in layer_iters]
        chunk_hasnext = np.ones(len(chunks), bool)

        pac_min = min([l.index.get_level_values("seq.pac")[0] for a,l in chunks])
        #seq_coords = get_full_coords(pac_min)


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

                pacs = chunk_layers.index.get_level_values("seq.pac").unique()
                while len(pacs) <= 1:
                    alns, layers = next(layer_iters[i], (None, None))
                    if alns is None:
                        chunk_hasnext[i] = False
                        break

                    chunk_alignments = pd.concat([chunks[i][0], alns])
                    chunk_alignments = chunk_alignments[~chunk_alignments.index.duplicated()]
                    chunk_layers = pd.concat([chunk_layers, layers])
                    chunks[i] = (chunk_alignments,chunk_layers)
                    pacs = chunk_layers.index.get_level_values("seq.pac").unique()

                idx = chunk_layers.index
                if all_fwd or all_rev:
                    fwds = idx.get_level_values("seq.fwd").unique()
                    all_fwd = all_fwd and np.all(fwds == 1)
                    all_rev = all_rev and np.all(fwds == 0)
                pacs = idx.get_level_values("seq.pac").unique()

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

                #seq_coords, coords = next_coords(seq_coords, chunk_pacs, fwd, chunk_hasnext[i])
                #coords.set_kmers(self.index.mposs_to_kmers(coords.mposs, self.conf.is_rna, False))

                #TODO probably don't need masks anymore
                #masks = [
                #    l.index.get_level_values("seq.pac").isin(coords.pacs) & (l.index.get_level_values("seq.fwd") == fwd)
                #    for a,l in chunks]


                #layers = pd.concat([l[mask] for (a,l),mask in zip(chunks, masks)])
                #alns = pd.concat([a.loc[l[mask].index.droplevel(["seq.fwd","seq.pac"]).unique()] for (a,l),mask in zip(chunks, masks)])

                layers = pd.concat({(i+1) : l for i,(a,l) in enumerate(chunks)}, names=["track.name","seq.fwd","seq.pac","aln.id"])
                alns = pd.concat({(i+1) : a.loc[l.index.droplevel(["seq.fwd","seq.pac"]).unique()] for i,(a,l) in enumerate(chunks)}, names=["track.name","aln.id"])
                aln_ids = alns.index

                refs = self.index.pac_to_pos(layers.index.get_level_values("seq.pac"))
                coords = RefCoord(alns.iloc[0]["ref_name"], refs.min(), refs.max()+1)

                ret = self._tables_to_tracks(coords, alns, layers)

                #leftover collection
                for i in range(len(chunks)):
                    chunk_alns, chunk_layers = chunks[i]

                    #if last_chunk:
                    chunk_alns = chunk_alns.iloc[:0]
                    chunk_layers = chunk_layers.iloc[:0]
                    #else:
                    #chunk_layers = chunk_layers[~masks[i]]
                    #if len(chunk_layers) > 0:
                    #    ids = chunk_layers.index.droplevel(["seq.fwd","seq.pac"]).unique()
                    #    chunk_alns = chunk_alns.loc[ids]
                    #else:
                    #    chunk_alns = chunk_alns.iloc[:0]

                    chunks[i] = (chunk_alns,chunk_layers)

                yield ret

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None, ignore_bam=False):

        if read_filter is None and self.read_index is not None:
            read_filter = self.read_index.read_filter
        
        if ref_bounds is not None and not isinstance(ref_bounds, RefCoord):
            ref_bounds = RefCoord(ref_bounds)

        if (self.coords is None or
            (read_filter is not None and 
             len(self.get_all_reads().intersection(read_filter)) < len(read_filter)) or
            (ref_bounds is not None and not self.coords.contains(ref_bounds))):
            gen = self.iter_reads_db(read_filter, ref_bounds, full_overlap, max_reads, ignore_bam)
        else:
            gen = self.iter_reads_slice(read_filter, ref_bounds)

        return gen
        #for read_id,chunk in gen:
        #    yield read_id,chunk
            
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

    def iter_reads_db(self, reads, ref_bounds, full_overlap, max_reads, ignore_bam=False):
        if ref_bounds is not None:
            self._set_ref_bounds(ref_bounds)
        if reads is None:
            reads = self.prms.read_filter
        if max_reads is None:
            max_reads = self.prms.max_reads

        #TODO handle multiple inputs properly
        for io in self.inputs:
            if ignore_bam and io == self.bam_in:
                continue

            aln_iter = io.iter_alns(
                self.db_layers, 
                self._aln_track_ids,
                self.coords,
                read_id=reads,
                full_overlap=full_overlap,
                ref_index=self.index)

            for aln in aln_iter:
                yield aln
        
    def _tables_to_tracks(self, coords, alignments, layers):
        tracks = dict()
        self.new_alignment = False
        self.new_layers = set()

        layer_alns = layers.index.get_level_values("aln.id")

        if self.prms.shared_refs_only or self.prms.min_coverage > 1:
            #track_covs = alignments.loc[layer_alns, ["track.name"]] 
            #track_covs = track_covs.set_index(layers.index)
            #track_covs = track_covs.reset_index("aln.id", drop=True) 
            #track_covs = track_covs.set_index("track.name", append=True) 
            track_covs = layers.index.droplevel("aln.id")
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

            l = layers.index.drop
            layers = layers[layers.index.droplevel(["track.name","aln.id"]).isin(idx)] # .loc[(slice(None),idx.get_level_values(0),idx.get_level_values(1),slice(None))]
            layer_alns = layers.index.droplevel(["seq.fwd","seq.pac"])
            alignments = alignments.loc[layer_alns.unique()]

        aln_groups = alignments.index.unique("track.name")
        for parent in self.alns:
            if parent.id in aln_groups:
                track_alns = alignments.loc[parent.id]
                track_layers = layers.loc[parent.id].droplevel("seq.fwd")
            else:
                track_alns = alignments.iloc[:0] 
                track_layers = layers.iloc[:0]   

            track = AlnTrack(parent, coords)#, track_alns, track_layers)

            #if not track.empty:
            #track.calc_layers(self.fn_layers)

            tracks[parent.name] = track

        tracks = Tracks(self, coords)
        tracks.layers = layers #, axis=0, names=["track.name", "aln.id", "seq.pos"])
        tracks.alignments = alignments #, axis=0, names=["track", "id"]).sort_index()#.sort_values(["fwd","ref_start"])

        #if not tracks.all_empty:
        #    tracks.load_compare(alignments.index.to_numpy())
        tracks.calc_refstats()

        return tracks

    def close(self):
        for io in self.inputs:
            io.close()
        if self.output is not None and (len(self.inputs) == 0 or self.output != self.inputs[0]):
            self.output.close()
