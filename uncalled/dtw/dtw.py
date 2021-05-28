#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict
import re
import time
from typing import NamedTuple
from matplotlib.colors import Normalize
import pandas as pd

from ..pafstats import parse_paf
from ..config import Config
from _uncalled import PORE_MODELS

class Track:
    CONF_FNAME = "conf.toml"
    ALN_DIR = "alns"
    ALN_SUFFIX = ".pkl"

    INDEX_FNAME = "filename_mapping.txt"
    INDEX_HEADER = "read_id\tfilename\taln_file"

    #REF_LOCS = "read_id\tsample_start\tsample_end\tref_name\tref_start\tref_end\tstrand"

    #TODO make fast5_file alias for filename in Fast5Reader
    #also, in fast5 reader make static fast5 filename parser

    WRITE_MODE = "w"
    READ_MODE = "r"
    MODES = {WRITE_MODE, READ_MODE}

    def __init__(self, path, mode="r", conf=None, overwrite=False):
        self.path = path.strip("/")
        self.mode = mode

        self.conf = conf if conf is not None else Config()

        if mode == self.WRITE_MODE:
            os.makedirs(self.aln_dir, exist_ok=overwrite)

        self.index_file = open(self.index_filename, mode)

        if mode == self.READ_MODE:
            self.conf.load_toml(self.config_filename)
            self.conf.fast5_reader.fast5_index = self.index_filename
            self._load_index()
            #self._load_reads()

        elif mode == self.WRITE_MODE:
            self.conf.to_toml(self.config_filename)
            self.index_file.write(self.INDEX_HEADER + "\n")

    @property
    def read_ids(self):
        return list(self.index.index)

    def add_read(self, read_id, fast5_fname, rae_df):
        if self.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")

        aln_fname = self.aln_fname(read_id)
        rae_df.to_pickle(aln_fname)

        self.index_file.write("\t".join([read_id, fast5_fname, aln_fname]) + "\n")

    def _load_index(self):
        self.index = pd.read_csv(self.index_file, sep="\t", index_col="read_id")

    def get_aln(self, read_id, ref_bounds=None):
        aln = pd.read_pickle(self.aln_fname(read_id))

        if ref_bounds is not None:
            _,st,en = ref_bounds[:3]
            return aln[st:en]
        return aln

    def get_matrix(self, ref_bounds, mm2_paf=None, partial_overlap=False):

        if mm2_paf is None:
            mm2_paf = self.conf.align.mm2_paf

        read_filter = set(self.index.index)
        mm2s = {p.qr_name : p
                 for p in parse_paf(
                    mm2_paf,
                    ref_bounds,
                    read_filter=read_filter,
                    full_overlap=not partial_overlap
        )}

        mat = TrackMatrix(self, ref_bounds, mm2s, conf=self.conf)

        for read_id,read in self.index.iterrows():
            mm2 = mm2s.get(read_id, None)
            if mm2 is None: continue
            df = pd.read_pickle(read["aln_file"])
            mat._add_read(df, mm2)

        mat._flatten()

        return mat

    def close(self):
        self.index_file.close()

    @property
    def config_filename(self):
        return os.path.join(self.path, self.CONF_FNAME)

    @property
    def index_filename(self):
        return os.path.join(self.path, self.INDEX_FNAME)
    
    @property
    def aln_dir(self):
        return os.path.join(self.path, self.ALN_DIR)
    
    def aln_fname(self, read_id):
        return os.path.join(self.aln_dir, read_id+self.ALN_SUFFIX)

class TrackMatrix:

    KMER_LAYER = 0
    PA_LAYER = 1
    DWELL_LAYER = 2
    PA_DIFF_LAYER = 3

    def __init__(self, track, ref_bounds, mm2s, height=None, conf=None):
        self.conf = conf if conf is not None else Config()
        self.track = track

        self.ref_bounds = ref_bounds
        self.width = self.ref_end-self.ref_start
        self.height = height

        self._layers = defaultdict(list)
        self.reads = defaultdict(list)
        self.mask = list()
        self.mm2s = dict()

        model_name = self.conf.mapper.pore_model
        #TODO probably need to rethink fwd/rev compl, but either way clean this up
        if model_name.endswith("_compl"):
            model_name = model_name[:-5]+"templ"

        self.model = PORE_MODELS[model_name]

        self.ref_to_x = pd.Series(
                np.arange(self.width, dtype=int),
                index=pd.RangeIndex(self.ref_start, self.ref_end)
        )

    def __getitem__(self, i):
        return self._layers.__getitem__(i)
    
    @property
    def kmer(self):
        return self._layers[self.KMER_LAYER]

    @property
    def pa(self):
        return self._layers[self.PA_LAYER]
    
    @property
    def dwell(self):
        return self._layers[self.DWELL_LAYER]
    
    @property
    def pa_diff(self):
        return self._layers[self.PA_DIFF_LAYER]

    def _add_read(self, df, mm2_paf):

        dtw_roi = df.loc[self.ref_start:self.ref_end-1]

        xs = self.ref_to_x.reindex(self.ref_to_x.index.intersection(dtw_roi.index))

        roi_mask = np.zeros(self.width)
        roi_mask[xs] = True
        self.mask.append(roi_mask)

        pa_diffs = self.model.match_diff(dtw_roi['mean'], dtw_roi['kmer'])
        dwell = 1000 * dtw_roi['length'] / self.conf.read_buffer.sample_rate

        self._add_layer_row(self.KMER_LAYER, dtw_roi['kmer'], xs)
        self._add_layer_row(self.PA_LAYER, dtw_roi['mean'], xs)
        self._add_layer_row(self.PA_DIFF_LAYER, pa_diffs, xs)
        self._add_layer_row(self.DWELL_LAYER, dwell, xs)

        self.reads['ref_start'].append(df.index.min())
        self.reads['id'].append(mm2_paf.qr_name)
        self.reads['fwd'].append(mm2_paf.is_fwd)

        self.mm2s[mm2_paf.qr_name] = mm2_paf

    def _add_layer_row(self, layer, vals, xs):
        row = np.zeros(self.width)
        row[xs] = vals
        self._layers[layer].append(row)
        return row

    def _flatten(self):
        self.reads = pd.DataFrame(self.reads) \
                     .sort_values(['fwd', 'ref_start'], ascending=[False, True])

        self.has_fwd = np.any(self.reads['fwd'])
        self.has_rev = not np.all(self.reads['fwd'])

        read_order = self.reads.index.to_numpy()

        layer_names = sorted(self._layers.keys())

        self._layers = np.stack([
            np.stack(self._layers[l])[read_order] for l in layer_names
        ])
        self.mask = np.stack(self.mask)[read_order].astype(bool)

        self.height = len(self.reads)

        self.norms = [Normalize(np.min(l[self.mask]), np.max(l[self.mask])) for l in self._layers]
        for l in [self.PA_LAYER, self.DWELL_LAYER]:
            lmask = self._layers[l][self.mask]
            self.norms[l].vmax = min(
                lmask.max(),
                np.median(lmask) + 2 * lmask.std()
            )

    def sort(self, layer, ref):
        order = np.argsort(self._layers[layer,:,ref])
        self._layers = self._layers[:,order,:]
        self.reads = self.reads.iloc[order]

    @property
    def ref_name(self):
        return self.ref_bounds[0]

    @property
    def ref_start(self):
        return self.ref_bounds[1]

    @property
    def ref_end(self):
        return self.ref_bounds[2]

def ref_coords(coord_str):
    spl = coord_str.split(":")
    ch = spl[0]
    st,en = spl[1].split("-")

    coord = (ch, int(st), int(en))

    if len(spl) == 2:
        return coord
    else:
        return coord + (spl[2] == "+",)


