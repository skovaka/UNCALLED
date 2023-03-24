import os
from collections import Sequence

import sys
import numpy as np
import pandas as pd
import h5py

from _uncalled import PoreModelParams, ArrayU32, ArrayU16
import _uncalled
from . import config

PORE_MODEL_PRESETS = {}

WORKFLOW_PRESETS = {
   "dna_r10.3_450bps" : "",
   "dna_r10.3_450bpsm" : "",
   "dna_r10.4.1_e8.2_260bps" : "",
   "dna_r10.4.1_e8.2_400bps" : "",
   "dna_r10_450bps" : "",
   "dna_r10.4_e8.1" : "",
   "dna_r10.4_e8.1m" : "",
   "dna_r9.4.1_450bps" : "r94_dna",
   "dna_r9.4.1_e8.1" : "r94_dna",
   "dna_r9.4.1_e8.1m" : "r94_dna",
   "dna_r9.5_450bps" : "r94_dna",
   "rna_r9.4.1_70bps" : "r94_rna"
}

CACHE = dict()

PARAM_TYPES = {
    "name" :  str,
    "k" : np.int32,
    "shift" :   np.int32,
    "pa_mean" : np.float32, 
    "pa_stdv" : np.float32, 
    "norm_max" : np.float32,
    "reverse" : bool,
    "complement" : bool
}

class PoreModel:

    @staticmethod
    def _param_defaults():
        return PoreModelParams(config._DEFAULTS.pore_model)

    @staticmethod
    def get_kmer_shift(k):
        return (k - 1) // 2

    #TODO load params like normal, maybe need python wrapper
    #store in pore model comments or binary
    #usually pass params, optionaly df/cache
    def __init__(self, *args, model=None, df=None, extra_cols=True, cache=True, **kwargs):
        self.conf, prms = config._init_group(
            "pore_model", _param_names=["name", "k", "shift", "norm_max", "reverse", "complement"], *args, **kwargs)

        self.instance = None

        is_preset = False

        self._cols = dict()

        vals = None

        self.extra_cols = extra_cols
        self.extra = pd.DataFrame() if self.extra_cols else None

        is_preset = False

        if model is not None: 
            if isinstance(getattr(model, "PRMS", None), PoreModelParams):
                if isinstance(model, PoreModel):
                    self._init(model.PRMS, model.instance)
                else:
                    self._init(model.PRMS, model)
            
            elif isinstance(model, pd.DataFrame):
                vals = self._vals_from_df(prms, model, True)
                self._init_new(prms, *vals)

            elif isinstance(model, dict):
                vals = self._vals_from_dict(prms, model)
                self._init_new(prms, *vals)

            else:
                raise TypeError(f"Invalid PoreModel type: {type(model)}")

        elif len(prms.name) > 0:

            cache_key = prms.to_tuple()

            if cache and cache_key in CACHE:
                self._init(prms, CACHE[cache_key])

            elif os.path.exists(prms.name):
                ext = prms.name.split(".")[-1]
                loader = self.FILE_LOADERS.get(ext, PoreModel._vals_from_tsv)
                vals = loader(self, prms)
                self._init_new(prms, *vals)
            
                if cache:
                    CACHE[cache_key] = self.instance

            else:
                models = ", ".join(PORE_MODEL_PRESETS.keys())
                raise FileNotFoundError(f"PoreModel file not found: {prms.name}")
        else:
            self._init_new(prms)

    def _init_new(self, prms, *args):
        if prms.k <= 8:
            ModelType = _uncalled.PoreModelU16
        else:
            ModelType = _uncalled.PoreModelU32

        if prms.shift < 0:
            prms.shift = PoreModel.get_kmer_shift(prms.k)


        self._init(prms, ModelType(prms, *args))

    def _init(self, prms, instance):
        self.instance = instance
        self.ModelType = type(instance)

        if prms.name is not None:
            self.PRMS.name = prms.name
        if prms.k > 0:
            self.PRMS.k = prms.k
        if prms.shift >= 0:
            self.PRMS.shift = prms.shift

        if self.K > 8:
            self.kmer_dtype = "uint32"
            self.array_type = ArrayU32
        else:
            self.kmer_dtype = "uint16"
            self.array_type = ArrayU16

        self.KMERS = np.arange(self.KMER_COUNT)
        #self.KMER_STRS = self.kmer_to_str(self.KMERS)
        self._KMER_STRS = None

        self._cols["current.mean"] = self.current.mean.to_numpy()
        self._cols["current.stdv"] = self.current.stdv.to_numpy()

    @property
    def KMER_STRS(self):
        if self._KMER_STRS is None:
            print("KEMRS")
            self._KMER_STRS = self.kmer_to_str(self.KMERS)
        return self._KMER_STRS

    @property
    def name(self):
        return self.PRMS.name

    @property
    def kmer_trim(self):
        return (self.PRMS.shift, self.K-self.PRMS.shift-1)

    COLUMNS = {"kmer", "current.mean", "current.stdv"}
    TSV_RENAME = {
        "current" : "current.mean",
        "mean" : "current.mean", 
        "stdv" : "current.stdv",
        "level_mean" : "current.mean", 
        "level_stdv" : "current.stdv",
        "sd"         : "current.stdv"
    }

    def _vals_from_df(self, prms, df, preprocess):
        df = df.rename(columns=self.TSV_RENAME)
        extra = df.columns.difference(self.COLUMNS)
        for col in extra:
            self._cols[col] = df[col].to_numpy()

        if preprocess:
            df = df.reset_index().sort_values("kmer")

        if prms.k < 0:
            kmer_lens = df["kmer"].str.len().value_counts()
            if len(kmer_lens) > 1:
                raise ValueError("All kmer lengths must be the same, found lengths: " + ", ".join(map(str, kmer_lens.index)))
            prms.k = kmer_lens.index[0]

        if "current.stdv" in df:
            stdv = df["current.stdv"].to_numpy()
        else:
            stdv = []
        return (df["current.mean"].to_numpy(), stdv, preprocess)
    
    def _usecol(self, name):
        return self.extra_cols or name in self.COLUMNS or name in self.TSV_RENAME

    def _vals_from_dict(self, prms, d):
        for name,typ in PARAM_TYPES.items():
            dname = "_"+name
            if not dname in d:
                sys.stderr.write(f"Warning: missing PoreModel parameter {name}\n")
            else:
                setattr(prms, name, typ(d[dname]))
                del d[dname]

        for k,v in d.items():
            if k in self.TSV_RENAME:
                k = self.TSV_RENAME[k]
            if k not in self.COLUMNS:
                self._cols[k] = v

        return (d["current.mean"], d["current.stdv"], False)

        #return self._vals_from_df(prms, pd.DataFrame(d), False)

    def _vals_from_npz(self, prms):
        d = dict(np.load(prms.name))
        return self._vals_from_dict(prms, d)

    def _vals_from_tsv(self, prms):
        df = pd.read_csv(prms.name, sep="\s+", comment="#", usecols=self._usecol)
        return self._vals_from_df(prms, df, True)

    def _vals_from_hdf5(self, prms):
        handle = h5py.File(prms.name, "r")
        df = pd.DataFrame(handle["model"][()]).reset_index()
        return self._vals_from_df(prms, df, True)

    def keys(self):
        return self._cols.keys()

    def __getitem__(self, idx):
        #TODO store master self._fields = {"mean" : self.means, ..., extra...}
        #no, maybe just wrap in series
        if isinstance(idx, str):
            return self._cols[idx]
        return self.mean.current[self.kmer_array(idx)]

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise ValueError(f"PoreModel has no attribute '{name}'")
        return self.instance.__getattribute__(name)

    def kmer_array(self, kmer):
        arr = np.array(kmer)
        if arr.shape == ():
            arr.shape = 1

        if arr.dtype.type in {np.str_, np.bytes_}:

            #TODO add option to fully check BP validity
            if not np.all(np.char.str_len(arr) == self.K):
                raise RuntimeError("All k-mers must be %d bases long" % self.K)

            arr = np.array([self.str_to_kmer(k) for k in arr])
        return self.array_type(arr.astype(self.kmer_dtype))
        #return arr

    def str_to_kmer(self, kmer):
        if isinstance(kmer, (Sequence, np.ndarray, pd.Series, self.array_type, list, tuple)):
            return np.array([self.instance.str_to_kmer(k, 0) for k in kmer])
        return self.instance.str_to_kmer(kmer, 0)
            

    def kmer_to_str(self, kmer, dtype=str):
        #, self.ModelType.KmerArray
        if isinstance(kmer, (Sequence, np.ndarray, pd.Series, self.array_type)):
            return self.instance.kmer_to_arr(kmer).astype(dtype)
        return dtype(self.instance.kmer_to_str(kmer))

    def norm_pdf(self, current, kmer):
        return self.instance.norm_pdf(self, current, self.kmer_array(kmer))

    def abs_diff(self, current, kmer):
        return self.instance.abs_diff(self, current, self.kmer_array(kmer))

    def to_dict(self, kmer_str=True, params=False):
        d = self._cols.copy()
        d.update(self.extra)
        if kmer_str:
            d["kmer"] = self.KMER_STRS
        else:
            d["kmer"] = self.KMERS
        if params:
            for name in PARAM_TYPES.keys():
                d[f"_{name}"] = getattr(self.PRMS, name)
        return d

    def to_df(self, kmer_str=True):
        return pd.DataFrame(self.to_dict(kmer_str))

    def to_tsv(self, out=None):
        return self.to_df().to_csv(out, sep="\t", index=False)

    def to_npz(self, fname):
        np.savez(fname, **self.to_dict(params=True))

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

        means = self.current.mean.to_numpy() * scale + shift
        vals = np.ravel(np.dstack([means,self.stdvs]))
        return PoreModel(self.ModelType(vals), name=self.name)
    
    
    def __setstate__(self, d):
        self.__init__(model=d)

    def __getstate__(self):
        return self.to_dict(params=True)
    
    def __repr__(self):
        ret = "<PoreModel mean=%.3f stdv=%.3f>\n" % (self.model_mean, self.model_stdv)
        ret += str(self.to_df())
        return ret[:-1]
    
    FILE_LOADERS = {
        "npz" : _vals_from_npz,
        "h5" : _vals_from_hdf5,
        "hdf5" : _vals_from_hdf5,
    }

