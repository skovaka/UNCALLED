#!/usr/bin/env python3

import sys, os
import sqlite3
import numpy as np
import pandas as pd
import collections
import time
from collections import defaultdict
import scipy

from .db import TrackSQL, delete, Eventalign, INPUT_PARAMS, OUTPUT_PARAMS
from .aln_track import AlnTrack, LAYERS, parse_layers
from ..index import load_index, RefCoord, str_to_coord
from ..pore_model import PoreModel
from ..fast5 import Fast5Reader, parse_read_ids
from .. import config


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
        if len(args) > 0 and isinstance(args[0], Tracks):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

    def _init_new(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("tracks", copy_conf=True, *args, **kwargs)

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers, add_deps=False))

        self.model = None#PoreModel(self.conf.pore_model) 

        #TODO refactor into TrackSQL, abstract into self.io

        self.alns = list()
        self.output_tracks = dict()
        self._tracks = dict()
        
        self._init_io()

        if self.prms.index_prefix is None:
            raise RuntimeError("Failed to load reference index")

        self._aln_track_ids = [t.id for t in self.alns]

        self.index = load_index(self.model.K, self.prms.index_prefix)

        self.coords = self._ref_bounds_to_coords(self.prms.ref_bounds)

        if self.coords is not None:
            self.load()

        if self.prms.load_fast5s and self.input is not None:
            fast5_reads = list()
            fast5_reads.append(self.input.get_fast5_index(self._aln_track_ids))
            fast5_reads = pd.concat(fast5_reads)
            files = fast5_reads["filename"].unique()
            self.fast5s = Fast5Reader(
                index=fast5_reads, 
                conf=self.conf)

            reads = fast5_reads["read_id"]
            self.reads = pd.Index(reads[~reads.duplicated()])

            for track in self.alns:
                track.fast5s = self.fast5s

    def _init_slice(self, parent, coords, tracks):
        self.conf = parent.conf 
        self.prms = parent.prms

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers, add_deps=False))

        self.output_tracks = parent.output_tracks
        self.index = parent.index
        self.fast5s = parent.fast5s
        self._aln_track_ids = parent._aln_track_ids
        self.refstats = None

        self.input = parent.input
        self.output = parent.output

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
            if LAYERS[group][layer].fn is None:
                if group in CMP_GROUPS:
                    self.cmp_layers.append((group, layer))
                else:
                    self.db_layers.append((group, layer))
            else:
                self.fn_layers.append((group, layer))

    def __len__(self):
        return len(self.alns)

    def __getitem__(self, i):
        return self.alns[i]

    def _init_io(self):
        in_prms = [getattr(self.prms.io, p) is not None for p in INPUT_PARAMS]
        out_prms = [getattr(self.prms.io, p) is not None for p in OUTPUT_PARAMS]
        if np.sum(in_prms) > 1:
            raise ValueError("No more than one input can be specified")
        if np.sum(out_prms) > 1:
            raise ValueError("No more than one output can be specified")

        if np.any(in_prms):
            in_format = INPUT_PARAMS[in_prms][0]
            if in_format == "db_in":
                self.input = TrackSQL(self.conf, self.model, "r")
                self.model = self.input.model
            elif in_format == "eventalign_in":
                raise ValueError("EVENTALGIN NOT SUPPORTED YET")
            elif in_format == "tombo_in":
                raise ValueError("TOMBO NOT SUPPORTED YET")

            for track in self.input.tracks:
                self._add_track(track.name, track)
        else:
            self.input = None

        if self.model is None:
            self.model = PoreModel(self.conf.pore_model.name)

        if np.any(out_prms):
            out_format = OUTPUT_PARAMS[out_prms][0]
            if out_format == "db_out":
                self.output = TrackSQL(self.conf, self.model, "w")
            elif out_format == "eventalign_out":
                self.output = Eventalign(self.conf, self.model, "w")

            for track in self.output.tracks:
                self.output_tracks[track.name] = track
        else:
            self.output = None

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
            return self.output.tracks[0]
            #if len(self.output_tracks) == 1:
            #    track_name, = self.output_tracks
            #else:
            #    raise ValueError("Must specify track name when using multiple output tracks")
        return self.output_tracks[track_name]

    def collapse_events(self, dtw):
        dtw["cuml_mean"] = dtw["length"] * dtw["current"]

        grp = dtw.groupby(level=0)

        lengths = grp["length"].sum()


        dtw = pd.DataFrame({
            "start"  : grp["start"].min().astype("uint32"),
            "length" : lengths.astype("uint32"),
            "current"   : grp["cuml_mean"].sum() / lengths,
            "kmer" : grp["kmer"].first()
        })

        return dtw.sort_index()

    def write_dtw_events(self, events=None, track_name=None, aln_id=None):
        if events is not None and self.output.FORMAT != "eventalign":
            events = self.collapse_events(events)
            overwrite = False
        else:
            overwrite = True

        self.write_layers("dtw", events, track_name, aln_id, overwrite)
        
    def write_layers(self, group, layers=None, track_name=None, aln_id=None, overwrite=False):
        track = self._track_or_default(track_name)

        if layers.index.names[0] == "mref":
            layers = layers.rename(index=self.index.mref_to_ref, level=0)
        elif layers.index.names[0] == "pac":
            layers = layers.rename(index=self.index.pac_to_ref, level=0)

        layers.index.names = ("ref",)

        if layers is not None:
            df = track.add_layer_group(group, layers, aln_id, overwrite)
        else:
            df = track.layers[group]

        self.output.write_layers(df.rename(index=track.coords.ref_to_pac, level=0))

    def write_alignment(self, read_id, fast5, coords, layers={}, track_name=None):
        track = self._track_or_default(track_name)

        aln_id = self.output.init_alignment(read_id, fast5)

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
        
        self.output.write_alignment(track.alignments)

        track.layers = None

        track.coords = coords

        for group,vals in layers.items():
            if group == "dtw":
                self.write_dtw_events(vals, track.name, aln_id)
            else:
                self.write_layers(group, vals, track.name, aln_id)

        return aln_id, coords

    @property
    def input_count(self):
        return len(self.alns)

    def _verify_read(self):
        if len(self._aln_track_ids) == 0:
            raise RuntimeError("No input tracks have been loaded")

    #TODO read_ids, track_alns, max_cov(?)
    def slice(self, ref_start=None, ref_end=None, reads=None, order=["fwd","ref_start"]):
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

        tracks = dict()
        for name,track in self._tracks.items():
            if isinstance(track, pd.DataFrame):
                tracks[name] = track.loc[coords.refs]

            elif isinstance(track, AlnTrack):
                tracks[name] = track.slice(coords, reads=reads, order=order)
                #if stranded:
                #    if tracks[name].all_fwd():



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

    def slice_shared_reads(self):
        return self.slice(reads=self.get_shared_reads(), order="read_id")

    def slice_full_overlap(self):
        read_ids = None
        rmin = self.coords.refs.min()
        rmax = self.coords.refs.max()
        for track in self.alns:
            alns = track.alignments[(track.alignments["ref_start"] <= rmin) & (track.alignments["ref_end"] >= rmax)]
            if read_ids is None:
                read_ids = pd.Index(alns["read_id"])
            else:
                read_ids = read_ids.intersection(alns["read_id"])

        return self.slice(reads=read_ids)
            
    def load(self, ref_bounds=None, full_overlap=None, read_filter=None, load_mat=False):
        self._verify_read()

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

        layers = self.input.query_layers(self.db_layers, self._aln_track_ids, self.coords, full_overlap=full_overlap, read_id=read_filter).droplevel(0)

        ids = layers.index.get_level_values("aln_id").unique().to_numpy()

        alignments = self.input.query_alignments(aln_id=ids)#self._aln_track_ids, coords=self.coords, full_overlap=full_overlap)


        for track in self.alns:
            track_alns = alignments[alignments["track_id"] == track.id]
            i = layers.index.get_level_values("aln_id").isin(track_alns.index)
            track_layers = layers.iloc[i]

            track.set_data(self.coords, track_alns, track_layers)
            track.calc_layers(self.fn_layers)
        
        self.load_compare(ids)
        self.calc_refstats()

        if load_mat:
            for track in self.alns:
                track.load_mat()

        return self.alns

    def load_compare(self, aln_ids=None):
        if len(self.cmp_layers) == 0:
            return

        self.cmp = self.input.query_compare(self.cmp_layers, self._aln_track_ids, self.coords, aln_ids)
        #TODO add aln_ids

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

            try:
                if bcalns is not None:
                    _add_group("bc_cmp", bcalns)
                if dtws is not None:
                    _add_group("cmp", dtws)
            except:
                sys.stderr.write("Failed to write compare group\n")
                sys.stderr.write(str(track.alignments))

    def calc_compare(self, group_b, calc_jaccard, calc_mean_ref_dist, save):
        if len(self.alns) > 0:
            alns = self.alns
        elif len(self.output_tracks) > 0:
            alns = list(self.output_tracks.values())
        else:
            raise ValueError("Must input at least one track")

        cols = alns[0].layers.columns.get_level_values("group").unique()
        if (group_b == "dtw" and "cmp" in cols) or (group_b == "bcaln" and "bc_cmp" in cols):
            sys.stderr.write(f"Read already has compare group. Skipping\n")
            return None

        if group_b == "dtw":
            if len(alns) != 2:
                raise ValueError("Must input exactly two tracks to compare dtw alignments")

            df = alns[0].cmp(alns[1], calc_jaccard, calc_mean_ref_dist)

        elif group_b == "bcaln":
            if len(alns) > 2:
                raise ValueError("Must input one or two tracks to compare dtw to bcaln")
            
            t = time.time()
            if len(alns) == 2:
                df = alns[0].bc_cmp(alns[1], calc_jaccard, calc_mean_ref_dist)
            else:
                df = alns[0].bc_cmp(None, calc_jaccard, calc_mean_ref_dist)

        df = df.dropna(how="all")

        if save:
            if "ref" in df.index.names:
                df.rename(index=lambda r: alns[0].coords.ref_to_pac(r), level=0, inplace=True)
            if self.output is None:
                output = self.input
            else:
                output = self.output
            output.write_layers( #TODO be explicit that output = input
                pd.concat({"cmp" : df}, names=["group", "layer"], axis=1), 
                index=["pac", "aln_a", "aln_b", "group_b"])
        #else:
        #    print(df.to_csv(sep="\t"))
        t = time.time()

    def calc_refstats(self, verbose_refs=False, cov=False):
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

        t0 = time.time()
        layer_iter = self.input.query_layers(
            self.db_layers, 
            self._aln_track_ids, 
            coords=coords, 
            order=["fwd","pac"],
            chunksize=self.prms.io.ref_chunksize)

        t0 = time.time()

        def get_full_coords(idx):
            fwd = idx[0]
            pac = idx[1]
            rid = self.index.pac_to_ref_id(pac)
            ref_len = self.index.get_ref_len(rid)
            ref_name = self.index.get_ref_name(rid)
            seq_refs = RefCoord(ref_name, 0, ref_len, fwd)
            return self.index.get_coord_space(seq_refs, self.conf.is_rna, load_kmers=False)

        def next_coords(seq_coords, idx):
            if len(idx) > 1:
                idx = idx[:-1]

            coords = seq_coords.pac_intersect(idx)
            if coords is None:
                seq_coords = get_full_coords(idx[0])
                coords = seq_coords.pac_intersect(idx)

            return seq_coords, coords

        chunk = next(layer_iter)
        seq_coords = get_full_coords(chunk.index[0])

        t0 = time.time()
                
        while len(chunk) > 0:

            if len(chunk.index.get_level_values("pac").unique()) <= 1:
                chunk = pd.concat([chunk, next(layer_iter, pd.DataFrame())])

            coord_idx = chunk.index.droplevel("aln_id").unique()
            #chunk_mrefs = chunk.index.get_level_values("pac").unique()

            seq_coords, coords = next_coords(seq_coords, coord_idx)

            
            coords.set_kmers(self.index.mrefs_to_kmers(coords.mrefs, self.conf.is_rna, False))

            mask = (chunk.index.get_level_values(0) == coords.fwd) & chunk.index.get_level_values(1).isin(coords.pacs)
            layers = chunk[mask]
            leftovers = chunk[~mask]

            aln_ids = layers.index.unique("aln_id").to_numpy()
            alns = self.input.query_alignments(self._aln_track_ids, aln_id=aln_ids)

            ret = self._tables_to_tracks(coords, alns, layers)

            chunk = leftovers

            if not ret.all_empty:
                yield ret

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None):
        
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

        t = time.time()
        
        layer_iter = self.input.query_layers(
            self.db_layers, 
            self._aln_track_ids,
            read_id=reads,
            coords=self.coords, 
            full_overlap=full_overlap, 
            order=["read_id", "fwd", "pac"],
            chunksize=self.prms.io.ref_chunksize)

        t = time.time()

        aln_leftovers = pd.DataFrame()
        layer_leftovers = pd.DataFrame()

        for layers in layer_iter:
            ids = layers.index \
                        .get_level_values("aln_id") \
                        .unique() \
                        .difference(aln_leftovers.index) \
                        .to_numpy()
            if len(ids) > 0:
                alignments = self.input.query_alignments(aln_id=ids)
            else:
                alignments = pd.DataFrame()

            alignments = pd.concat([aln_leftovers, alignments])
            layers = pd.concat([layer_leftovers, layers])

            aln_end = alignments["read_id"] == alignments["read_id"].iloc[-1]
            aln_leftovers = alignments.loc[aln_end]
            alignments = alignments.loc[~aln_end]

            layer_end = layers.index.get_level_values("aln_id").isin(aln_leftovers.index)
            layer_leftovers = layers.loc[layer_end]
            layers = layers.loc[~layer_end]
            layer_alns = layers.index.get_level_values("aln_id")

            for ref_name,ref_alns in alignments.groupby("ref_name"):
                coords = self._alns_to_coords(ref_alns)
                cache = self._tables_to_tracks(coords, ref_alns, layers)
                for ret in cache.iter_reads_slice():
                    yield ret

        if len(aln_leftovers) > 0:
            for ref_name,ref_alns in aln_leftovers.groupby("ref_name"):
                coords = self._alns_to_coords(ref_alns)
                cache = self._tables_to_tracks(coords, ref_alns, layer_leftovers)
                for ret in cache.iter_reads_slice():
                    yield ret

    def _tables_to_tracks(self, coords, alignments, layers):
        tracks = dict()

        layer_alns = layers.index.get_level_values("aln_id")

        if self.prms.shared_refs_only or self.prms.min_coverage > 1:
            track_covs = alignments.loc[layer_alns, ["track_id"]] \
                                   .set_index(layers.index) \
                                   .reset_index("aln_id", drop=True) \
                                   .set_index("track_id", append=True) \
                                   .index.value_counts()

            mask = track_covs >= self.prms.min_coverage
            if not np.any(mask):
                idx = pd.Index([])
            elif self.prms.shared_refs_only:
                track_counts = pd.MultiIndex.from_tuples(track_covs[mask].index) \
                                   .droplevel(2) \
                                   .value_counts()

                idx = pd.MultiIndex.from_tuples(track_counts.index[track_counts == len(self.alns)])
            else:
                idx = track_covs[mask].index.droplevel("aln_id").unique()

            layers = layers.loc[(idx.get_level_values(0),idx.get_level_values(1),slice(None)), slice(None)]
            layer_alns = layers.index.get_level_values("aln_id")
            alignments = alignments.loc[layer_alns.unique()]

        layers = layers.droplevel(0)

        aln_groups = alignments.groupby("track_id").groups
        for parent in self.alns:
            if parent.id in aln_groups:
                track_alns = alignments.loc[aln_groups[parent.id]]
                track_layers = layers.loc[layer_alns.isin(track_alns.index)]
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
        if self.input is not None:
            self.input.close()
        if self.output is not None:
            self.output.close()
