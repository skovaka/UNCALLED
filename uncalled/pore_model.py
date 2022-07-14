import os
from collections import Sequence

import sys
import numpy as np
import pandas as pd
import h5py

from _uncalled import PoreModelParams, ArrayU32, ArrayU16, PORE_MODEL_PRESETS
import _uncalled
from . import config

CACHE = dict()

class PoreModel:

    @staticmethod
    def _param_defaults():
        return PoreModelParams(config._DEFAULTS.pore_model)

    @staticmethod
    def get_kmer_shift(k):
        return (k - 1) // 2

    def __init__(self, model=None, name=None, reverse=None, complement=None, df=None, extra_cols=False, cache=True):
        is_preset = False

        self._cols = dict()

        if model is not None: 
            if isinstance(getattr(model, "PRMS", None), PoreModelParams):
                if isinstance(model, PoreModel):
                    self._init(model.prms, model.instance)
                else:
                    self._init(model.PRMS, model)
                return 

            if isinstance(model, str):

                if model in PORE_MODEL_PRESETS:
                    prms = PORE_MODEL_PRESETS[model].prms
                    is_preset = True
                else:
                    prms = self._param_defaults()
                    prms.name = model

            elif isinstance(model, PoreModelParams):
                prms = PoreModelParams(model)
                is_preset = model.name in PORE_MODEL_PRESETS

            else:
                raise TypeError("PoreModel model must be of type str, PoreModel, or PoreModel.Params")
        else:
            prms = self._param_defaults()
            if df is None: is_preset = True

        if reverse is not None: prms.reverse = reverse
        if complement is not None: prms.complement = complement

        self.extra_cols = extra_cols
        self.extra = pd.DataFrame() if self.extra_cols else None

        vals = None

        if df is not None:
            vals = self._vals_from_df(df)

        elif cache and prms.name in CACHE:
            sys.stderr.write(f"Using cached model: {prms.name}\n")
            self._init(prms, CACHE[prms.name])
            return

        elif not is_preset:
            if os.path.exists(prms.name):
                vals = self._vals_from_tsv(prms.name)
                #try:
                #except:
                #    try:
                #        vals = self._vals_from_hdf5(prms.name)
                #    except:
                #        raise ValueError("Unrecognized PoreModel file format. Must be a valid TSV or HDF5 file.")
            else:
                models = ", ".join(PORE_MODEL_PRESETS.keys())
                raise ValueError(
                    f"Unknown PoreModel: {prms.name}. Must be a filename, or one of {models}" \
                    
                )

        if vals is not None:
            prms.k = int(np.log2(len(vals)) / 2)
            self._init_new(prms, vals, prms.reverse, prms.complement)
        else:
            self._init_new(prms, prms)
            
        if cache:
            CACHE[prms.name] = self.instance

    def _init_new(self, prms, *args):
        ModelType = getattr(_uncalled, f"PoreModelK{prms.k}", None)
        if ModelType is None:
            raise ValueError(f"Invalid k-mer length {prms.k}")
        self._init(prms, ModelType(*args))

    def _init(self, prms, instance):
        self.instance = instance
        self.ModelType = type(instance)

        if prms.name is not None:
            self.PRMS.name = prms.name

        if self.K > 8:
            self.kmer_dtype = "uint32"
            self.array_type = ArrayU32
        else:
            self.kmer_dtype = "uint16"
            self.array_type = ArrayU16

        self.KMERS = np.arange(self.KMER_COUNT)
        self.KMER_STRS = self.kmer_to_str(self.KMERS)

        self._cols["mean"] = self.means
        self._cols["stdv"] = self.stdvs

    @property
    def name(self):
        return self.PRMS.name

    COLUMNS = {"kmer", "mean", "stdv"}
    TSV_RENAME = {
        "level_mean" : "mean", 
        "level_stdv" : "stdv",
        "sd"         : "stdv"
    }

    def _vals_from_df(self, df):
        df = df.rename(columns=self.TSV_RENAME)
        extra = df.columns.difference(self.COLUMNS)
        for col in extra:
            self._cols[col] = df[col].to_numpy()
        return np.ravel(df.sort_values("kmer")[["mean","stdv"]])
               

    def _usecol(self, name):
        return self.extra_cols or name in self.COLUMNS or name in self.TSV_RENAME

    def _vals_from_tsv(self, filename):
        df = pd.read_csv(filename, sep="\t", comment="#", usecols=self._usecol)
        return self._vals_from_df(df)

    def _vals_from_hdf5(self, filename):
        handle = h5py.File(filename, "r")
        df = pd.DataFrame(handle["model"][()]).reset_index()
        return self._vals_from_df(df)

    def keys(self):
        return self._cols.keys()

    def __getitem__(self, idx):
        #TODO store master self._fields = {"mean" : self.means, ..., extra...}
        #no, maybe just wrap in series
        if isinstance(idx, str):
            return self._cols[idx]
        return self.means[self.kmer_array(idx)]

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def kmer_array(self, kmer):
        arr = np.array(kmer)
        if arr.shape == ():
            arr.shape = 1

        if arr.dtype.type in {np.str_, np.bytes_}:

            #TODO add option to fully check BP validity
            if not np.all(np.char.str_len(arr) == self.K):
                raise RuntimeError("All k-mers must be %d bases long" % self.K)

            arr = self.str_to_kmer(arr)
        return self.array_type(arr.astype(self.kmer_dtype))
        #return arr

    def kmer_to_str(self, kmer, dtype=str):
        #, self.ModelType.KmerArray
        if isinstance(kmer, (Sequence, np.ndarray, pd.Series, self.array_type)):
            return self.ModelType.kmer_to_arr(kmer).astype(dtype)
        else:
            print(type(kmer), "NO")
        return dtype(self.ModelType.kmer_to_str(kmer))

    def norm_pdf(self, current, kmer):
        return self.instance.norm_pdf(self, current, self.kmer_array(kmer))

    def abs_diff(self, current, kmer):
        return self.instance.abs_diff(self, current, self.kmer_array(kmer))

    def to_df(self, kmer_str=True):
        df = pd.DataFrame(self._cols)

        if kmer_str:
            df["kmer"] = self.KMER_STRS 

        return df

    def to_tsv(self, out=None):
        return self.to_df().to_csv(out, sep="\t", index=False)

    def norm_mom_params(self, current, tgt_mean=None, tgt_stdv=None):
        tgt_mean = self.model_mean if tgt_mean is None else tgt_mean
        tgt_stdv = self.model_stdv if tgt_stdv is None else tgt_stdv
        scale = tgt_stdv / np.std(current)
        shift = tgt_mean - scale * np.mean(current)
        return scale, shift

    def norm_mad_params(self, current):
        shift = np.median(current)
        scale = 1 / np.median(np.abs(current - shift))
        shift *= -scale
        return scale, shift

    def get_normalized(self, scale, shift):
        #if len(args) == 2:
        #    tgt_mean,tgt_stdv = args
        #elif len(args) == 1:
        #    model = args[0]
        #    tgt_mean = model.model_mean
        #    tgt_stdv = model.model_stdv
        #scale,shift = self.norm_mom_params(self.means, tgt_mean, tgt_stdv)

        means = self.means * scale + shift
        vals = np.ravel(np.dstack([means,self.stdvs]))
        return PoreModel(self.ModelType(vals), name=self.name)
    
    
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

