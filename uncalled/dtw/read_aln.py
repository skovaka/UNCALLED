#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict, namedtuple
import re
import time
from matplotlib.colors import Normalize
import pandas as pd
import scipy.stats
import copy

import matplotlib.pyplot as plt

from ..pafstats import parse_paf, PafEntry
from ..config import Config
from ..argparse import Opt
from .. import nt, PoreModel


class ReadAln:

    def __init__(self, track_id, read_id, mrefs=None, index=None, is_rna=False):
        self.id = None
        self.track_id = track_id
        self.read_id = read_id
        self.index = index
        self.is_rna = is_rna

        self.dfs = set()

        self.mrefs = mrefs
        if mrefs is not None:
            self.set_coords(mrefs)

    def set_coords(self, coords):
        loc = self.index.mref_to_ref_bound(self.mref_start,self.mref_end,not self.is_rna)
        self.ref_bounds = RefCoord(loc.ref_name, loc.start, loc.end, loc.fwd)
        self.ref_id = loc.ref_id

    @property
    def empty(self):
        return not hasattr(self, "dtw") or len(self.dtw) == 0

    @property
    def mref_start(self):
        return self.mrefs.min()

    @property
    def mref_end(self):
        return self.mrefs.max()+1
    
    @property
    def ref_start(self):
        return self.ref_bounds.start

    @property
    def ref_end(self):
        return self.ref_bounds.end

    @property
    def ref_name(self):
        return self.ref_bounds.name

    #TODO ReadAln stores RefCoords, handles all this conversion?

    def mref_to_ref(self, mref):
        return self.index.mref_to_ref(mref) 

    def ref_to_mref(self, ref):
        return self.index.ref_to_mref(self.ref_name, ref, ref, self.is_fwd, self.is_rna)[0]

    def ref_to_samp(self, ref):
        return self.mref_to_samp(self.ref_to_mref(ref))
        
    def mref_to_samp(self, mref):
        i = np.clip(self.dtw['mref'].searchsorted(mref), 0, len(self.dtw)-1)
        return self.dtw['sample'][i]
    
    def calc_mref(self):
        self.dtw["mref"] = self.index.ref_to_mref(self.ref_id, self.dtw.index, self.is_fwd, self.is_rna)

    @property
    def is_fwd(self):
        return self.ref_bounds.fwd

    def sort_mref(self):
        self.dtw.sort_values("mref", inplace=True)

    def get_samp_bounds(self):
        samp_min = int(self.dtw['start'].min())
        max_i = self.dtw['start'].argmax()
        samp_max = int(self.dtw['start'].iloc[max_i] + self.dtw['length'].iloc[max_i])
        return samp_min, samp_max
    
    #def set_bands(self, bands):
    def set_df(self, df, name):
        if self.mrefs is None:
            self.mrefs = df.index
            self.set_coords(df.index)
        else:
            index = df.index.intersection(self.mrefs)
            df = df.reindex(index=index, copy=False)

        self.dfs.add(name)
        if not "aln_id" in df:
            df["aln_id"] = self.id

        setattr(self, name, df)

    def set_dtw(self, df):
        self.set_df(df, "dtw")

    def set_bcerr(self, df):
        self.set_df(df, "bcerr")

    def set_subevent_aln(self, dtw, ref_mirrored=False, kmer_str=False, ref_col="mref", start_col="start", length_col="length", mean_col="current", kmer_col="kmer"):

        dtw["cuml_mean"] = dtw[length_col] * dtw[mean_col]

        grp = dtw.groupby(ref_col)

        if kmer_col in dtw:
            if kmer_str:
                kmers = [nt.kmer_rev(nt.str_to_kmer(k,0)) for k in grp[kmer_col].first()]
            else:
                kmers = grp[kmer_col].first()
        else:
            kmers = None

        if ref_col == "mref":
            mrefs = grp[ref_col].first()
        else:
            mrefs = self.ref_to_mref(grp[ref_col].first())

        lengths = grp[length_col].sum()

        dtw = pd.DataFrame({
            "mref"    : mrefs.astype("int64"),
            "start"  : grp[start_col].min().astype("uint32"),
            "length" : lengths.astype("uint32"),
            "current"   : grp["cuml_mean"].sum() / lengths
        })

        if kmers is not None:
            dtw["kmer"] = kmers.astype("uint16")

        dtw = dtw.set_index("mref").sort_index()
        self.set_dtw(dtw)

    def get_index_kmers(self, index, kmer_shift=4):
        """Returns the k-mer sequence at the alignment reference coordinates"""
        #start = max(0, self.mref_start - kmer_shift)
        start = self.mref_start

        return pd.Series(
            index.get_kmers(self.mref_start, self.mref_end, self.is_rna),
            pd.RangeIndex(self.mref_start+nt.K-1, self.mref_end)
        )
        #kmers = index.get_kmers(start, self.mref_end, self.is_rna)
        #return np.array(kmers)

