import sys
import time
import glob
import os
import re
from collections import Sequence

import numpy as np
import pandas as pd
import h5py

from _uncalled import PoreModelK5, PoreModelK10, PoreModelParams
from . import config

#class PoreModel(_PoreModel):
class PoreModel:

    @staticmethod
    def _param_defaults():
        return PoreModelParams(config._DEFAULTS.pore_model)

    def __init__(self, model=None, name=None, reverse=None, complement=None, df=None):
        self.ModelType = PoreModelK5

        if model is not None: 
            if isinstance(model, PoreModelK5):
                self._init(name, model)
                return 

            if isinstance(model, PoreModelK10):
                self._init(name, model)
                return 

            if isinstance(model, str):
                prms = self._param_defaults()
                prms.name = model

            elif isinstance(model, PoreModelParams):
                prms = PoreModelParams(model)

            else:
                raise TypeError("PoreModel model must be of type str, PoreModel, or PoreModel.Params")
        else:
            prms = self._param_defaults()

        if reverse is not None: prms.reverse = reverse
        if complement is not None: prms.complement = complement

        #pd.DataFrame overrides
        if df is not None:
            vals = self._vals_from_df(df)

        elif self._init_preset(prms):
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
            models = ", ".join(PoreModelK5.get_preset_names())
            raise ValueError(
                f"Unknown PoreModel: {prms.name}. Must be a filename, or one of {models}" \
                
            )

        k = int(np.log2(len(vals)) / 2)
        if k == 5:
            self.ModelType = PoreModelK5
        elif k == 10:
            self.ModelType = PoreModelK10
        else:
            raise ValueError(f"Invalid k-mer length: {k}\n")

        self._init(name, vals, prms.reverse, prms.complement)

    def _init(self, name, *args):
        self.instance = self.ModelType(*args)
        if name is not None:
            self.PRMS.name = name

        self.KMERS = np.arange(self.KMER_COUNT)
        self.KMER_STRS = self.kmer_to_str(self.KMERS)


    def _init_preset(self, prms):
        if PoreModelK5.is_preset(prms.name):
            self.ModelType = PoreModelK5
        elif PoreModelK10.is_preset(prms.name):
            self.ModelType = PoreModelK10
        else:
            return False
        self._init(prms.name, prms)
        return True

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
        return self.means[self.kmer_array(kmer)]

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def kmer_array(self, kmer):
        arr = np.array(kmer)
        if arr.shape == ():
            arr.shape = 1

        if arr.dtype.type in {np.str_, np.bytes_}:

            #TODO add option to fully check BP validity
            if not np.all(np.char.str_len(arr) == K):
                raise RuntimeError("All k-mers must be %d bases long" % K)

            arr = str_to_kmer(arr)
        return self.KmerArray(arr.astype("uint16"))
        #return arr

    def kmer_to_str(self, kmer, dtype=str):
        if isinstance(kmer, (Sequence, self.ModelType.KmerArray, np.ndarray, pd.Series)):
            return self.ModelType.kmer_to_arr(kmer).astype(dtype)
        return dtype(self.ModelType.kmer_to_str(kmer))

    def norm_pdf(self, current, kmer):
        return self.instance.norm_pdf(self, current, self.kmer_array(kmer))

    def abs_diff(self, current, kmer):
        return self.instance.abs_diff(self, current, self.kmer_array(kmer))

    def to_df(self):
        return pd.DataFrame({
            "kmer" : self.KMER_STRS, 
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

        return(PoreModel(self.ModelType(vals), name=self.name), scale, shift)
    
    
    def __repr__(self):
        ret = "<PoreModel mean=%.3f stdv=%.3f>\n" % (self.model_mean, self.model_stdv)
        ret += "kmer    mean    stdv\n"
        def kmer_str(i):
            return "%s%8.3f%8.3f\n" % (self.KMER_STRS[i], self.means[i], self.stdvs[i])
        for i in self.KMERS[:3]:
            ret += kmer_str(i)
        ret += "...\n"
        for i in self.KMERS[-3:]:
            ret += kmer_str(i)
        return ret[:-1]

