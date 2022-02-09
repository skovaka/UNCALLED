import sys
import time
import glob
import os
import re
from collections import Sequence

import numpy as np
import pandas as pd
import h5py

from _uncalled import _PoreModel
from . import nt, config

class PoreModel(_PoreModel):

    @staticmethod
    def _param_defaults():
        return _PoreModel.Params(config._DEFAULTS.pore_model)

    def __init__(self, model=None, name=None, reverse=None, complement=None, df=None):
        if model is not None: 
            if isinstance(model, _PoreModel):
                self._init(name, model)
                return 

            if isinstance(model, str):
                prms = self._param_defaults()
                prms.name = model

            elif isinstance(model, PoreModel.Params):
                prms = _PoreModel.Params(model)

            else:
                raise TypeError("PoreModel model must be of type str, PoreModel, or PoreModel.Params")
        else:
            prms = self._param_defaults()

        if reverse is not None: prms.reverse = reverse
        if complement is not None: prms.complement = complement

        #pd.DataFrame overrides
        if df is not None:
            vals = _vals_from_df(df)

        elif self.is_preset(prms.name):
            self._init(name, prms)
            return

        elif os.path.exists(prms.name):
            try:
                vals = self._vals_from_tsv(prms.name)
            except:
                try:
                    vals = self._vals_from_hdf5(prms.name)
                except:
                    raise ValueError("Unrecognized PoreModel file format. Must be a valid TSV or HDF5 file.")
        else:
            raise ValueError(
                "PoreModel name must be a filename or one of {%s}" % \
                ", ".join(_PoreModel.get_preset_names())
            )

        self._init(name, vals, prms.reverse, prms.complement)

    def _init(self, name, *args):
        _PoreModel.__init__(self, *args)
        if name is not None:
            self.PRMS.name = name

    @property
    def name(self):
        return self.PRMS.name

    def _vals_from_df(self, df):
        return np.ravel(df.rename(columns=self.TSV_RENAME) \
                          .sort_values("kmer")[["mean","stdv"]])
               
    COLUMNS = {"kmer", "mean", "stdv"}
    TSV_RENAME = {
        "level_mean" : "mean", 
        "level_stdv" : "stdv",
        "sd"         : "stdv"
    }

    def _usecol(self, name):
        return name in self.COLUMNS or name in self.TSV_RENAME

    def _vals_from_tsv(self, filename):
        df = pd.read_csv(filename, sep="\t", comment="#", usecols=self._usecol)
        return self._vals_from_df(df)

    def _vals_from_hdf5(self, filename):
        handle = h5py.File(filename, "r")
        df = pd.DataFrame(handle["model"][()]).reset_index()
        return self._vals_from_df(df)

    def __getitem__(self, kmer):
        return self.means[nt.kmer_array(kmer)]

    def norm_pdf(self, current, kmer):
        return _PoreModel.norm_pdf(self, current, nt.kmer_array(kmer))

    def abs_diff(self, current, kmer):
        return _PoreModel.abs_diff(self, current, nt.kmer_array(kmer))

    def to_df(self):
        return pd.DataFrame({
            "kmer" : nt.KMER_STRS, 
            "mean" : self.means, 
            "stdv" : self.stdvs
        })

    def to_tsv(self, out=None):
        return self.to_df().to_csv(out, sep="\t", index=False)

    def norm_mom_params(self, current, tgt_mean=None, tgt_stdv=None):
        tgt_mean = self.model_mean if tgt_mean is None else tgt_mean
        tgt_stdv = self.model_stdv if tgt_stdv is None else tgt_stdv
        scale = tgt_stdv / np.std(current)
        shift = tgt_mean - scale * np.mean(current)
        return scale, shift

    def get_normalized(self, *args):
        if len(args) == 2:
            tgt_mean,tgt_stdv = args
        elif len(args) == 1:
            model = args[0]
            tgt_mean = model.model_mean
            tgt_stdv = model.model_stdv

        scale,shift = self.norm_mom_params(self.means, tgt_mean, tgt_stdv)
        new_means = self.means * scale + shift
        vals = np.ravel(np.dstack([new_means,self.stdvs]))

        return(PoreModel(_PoreModel(vals), name=self.name), scale, shift)
    
    
    def __repr__(self):
        ret = "<PoreModel mean=%.3f stdv=%.3f>\n" % (self.model_mean, self.model_stdv)
        ret += "kmer    mean    stdv\n"
        def kmer_str(i):
            return "%s%8.3f%8.3f\n" % (nt.KMER_STRS[i], self.means[i], self.stdvs[i])
        for i in nt.KMERS[:3]:
            ret += kmer_str(i)
        ret += "...\n"
        for i in nt.KMERS[-3:]:
            ret += kmer_str(i)
        return ret[:-1]

