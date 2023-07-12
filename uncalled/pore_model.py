import os
from collections import Sequence

import sys
import numpy as np
import pandas as pd
import h5py
import itertools

from _uncalled import PoreModelParams, ArrayU32, ArrayU16
import _uncalled
from . import config

PORE_MODEL_PRESETS = {}

CACHE = dict()

PARAM_TYPES = {
    "name" :  str,
    "k" : np.int32,
    "shift" :   np.int32,
    "pa_mean" : np.float32, 
    "pa_stdv" : np.float32, 
    "norm_max" : np.float32,
    "sample_rate" : np.float32,
    "bases_per_sec" : np.float32,
    "reverse" : bool,
    "complement" : bool,
    "flowcell" : str,
    "kit" : str
}


ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

class PoreModel:

    PRESET_DIR = os.path.join(ROOT_DIR, "models")
    PRESET_EXT = ".npz"

    PRESET_MAP = None
    PRESETS = {"dna_r10.4.1_400bps_9mer", "dna_r9.4_400bps_5mer", "rna_r9.4_70bps_5mer", "tombo/rna_r9.4_70bps_5mer"}

    @classmethod
    def _init_presets(cls):
        if cls.PRESET_MAP is None:
            df = pd.read_csv(
                os.path.join(cls.PRESET_DIR, "presets.tsv"), 
                sep="\t", index_col=("flowcell","kit")
            ).sort_index()

            cls.PRESET_MAP = df[df["preset_model"] != "_"]
        #if True:# cls.PRESETS is None:
        #    #cls.PRESETS = set()
        #    for root, dirs, files in os.walk(cls.PRESET_DIR):
        #        for fname in files:
        #            if fname.endswith(".npz"):
                        #cls.PRESETS.add(fname)

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

        self._base = dict()
        self._extra = dict() if extra_cols else None

        vals = None

        is_preset = False

        if model is not None: 
            if isinstance(getattr(model, "PRMS", None), PoreModelParams):
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
            cache_key = prms.to_key()
            if cache and cache_key in CACHE:
                self._init(prms, CACHE[cache_key])
                return

            if prms.name in self.PRESETS:
                filename = os.path.join(self.PRESET_DIR, prms.name + self.PRESET_EXT)
                ext = self.PRESET_EXT[1:]
            else:
                filename = prms.name
                ext = filename.split(".")[-1]

            if os.path.exists(filename):
                loader = self.FILE_LOADERS.get(ext, PoreModel._vals_from_tsv)
                vals = loader(self, filename, prms)
                self._init_new(prms, *vals)

            else:
                models = ", ".join(PORE_MODEL_PRESETS.keys())
                raise FileNotFoundError(f"PoreModel file not found: {filename}. Choose one of: {models}")

        else:
            self._init_new(prms)
            
        cache_key = self.PRMS.to_key()
        if cache and not cache_key in CACHE:
            CACHE[cache_key] = self

    def _init_new(self, prms, *args):
        if prms.k <= 8:
            ModelType = _uncalled.PoreModelU16
        else:
            ModelType = _uncalled.PoreModelU32

        if prms.shift < 0:
            prms.shift = PoreModel.get_kmer_shift(prms.k)

        self._init(prms, ModelType(prms, *args))

    def _init(self, prms, model):

        if isinstance(model, PoreModel):
            self.instance = model.instance
            self._base.update(model._base)
            self._extra.update(model._extra)
        else:
            self.instance = model
            self._base["current.mean"] = model.current.mean.to_numpy()
            self._base["current.stdv"] = model.current.stdv.to_numpy()
            self._base["current_sd.mean"] = model.current_sd.mean.to_numpy()
            self._base["current_sd.stdv"] = model.current_sd.stdv.to_numpy()

        self.ModelType = type(self.instance)

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
        self._KMER_STRS = None

    @property
    def KMER_STRS(self):
        if self._KMER_STRS is None:
            self._KMER_STRS = self.kmer_to_str(self.KMERS)
        return self._KMER_STRS

    @property
    def name(self):
        return self.PRMS.name

    @property
    def kmer_trim(self):
        return (self.PRMS.shift, self.K-self.PRMS.shift-1)

    @property
    def reverse(self):
        return self.PRMS.reverse

    COLUMNS = {"kmer", "current.mean", "current.stdv", "current_sd.mean", "current_sd.stdv"}
    TSV_RENAME = {
        "current" : "current.mean",
        "mean" : "current.mean", 
        "stdv" : "current.stdv",
        "level_mean" : "current.mean", 
        "level_stdv" : "current.stdv",
        "sd_mean"         : "current_sd.mean",
        "sd_stdv"         : "current_sd.stdv",
        "sd"         : "current.stdv",
    }

    def _vals_from_df(self, prms, df, preprocess):
        df = df.rename(columns=self.TSV_RENAME)
        if self._extra is not None:
            extra = df.columns.difference(self.COLUMNS)
            for col in extra:
                #self._cols[col] = df[col].to_numpy()
                self._extra[col] = df[col].to_numpy()

        if preprocess:
            df = df.reset_index().sort_values("kmer")

        if prms.k < 0:
            print(df["kmer"])
            kmer_lens = df["kmer"].str.len().value_counts()
            if len(kmer_lens) > 1:
                raise ValueError("All kmer lengths must be the same, found lengths: " + ", ".join(map(str, kmer_lens.index)))
            prms.k = kmer_lens.index[0]

        get = lambda c: df[c] if c in df else []

        return (get("current.mean"), get("current.stdv"), get("current_sd.mean"), get("current_sd.stdv"), preprocess)
    
    def _usecol(self, name):
        return self._extra is not None or name in self.COLUMNS or name in self.TSV_RENAME

    def _vals_from_dict(self, prms, d):
        for name,typ in PARAM_TYPES.items():
            dname = "_"+name
            if dname in d:
                old_val = getattr(prms, name)
                new_val = typ(d[dname])
                if ((typ != np.int32 or old_val < 0) and 
                    (not hasattr(new_val, "__len__") or len(new_val) > 0)):
                    setattr(prms, name, new_val)
                del d[dname]

        for k,v in d.items():
            if k in self.TSV_RENAME:
                k = self.TSV_RENAME[k]
            if k not in self.COLUMNS:
                self._base[k] = v

        get = lambda c: d[c] if c in d else []

        #return (d["current.mean"], get("current.stdv"], False)
        return (get("current.mean"), get("current.stdv"), get("current_sd.mean"), get("current_sd.stdv"), False)

        #return self._vals_from_df(prms, pd.DataFrame(d), False)

    def _vals_from_npz(self, filename, prms):
        d = dict(np.load(filename))
        return self._vals_from_dict(prms, d)

    def _vals_from_tsv(self, filename, prms):
        df = pd.read_csv(filename, sep="\s+", comment="#", usecols=self._usecol)
        return self._vals_from_df(prms, df, True)

    def _vals_from_hdf5(self, filename, prms):
        handle = h5py.File(filename, "r")
        df = pd.DataFrame(handle["model"][()])#.reset_index()
        return self._vals_from_df(prms, df, True)

    def keys(self):
        return itertools.chain(self._base.keys(), self._extra.keys())

    def __getitem__(self, idx):
        if isinstance(idx, str):
            return self._base[idx]
        return self.current.mean[self.kmer_array(idx)]

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            if hasattr(self.PRMS, name):
                return getattr(self.PRMS, name)
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

    #def norm_pdf(self, current, kmer):
    #    return self.instance.norm_pdf(current, self.kmer_array(kmer))

    def abs_diff(self, current, kmer):
        return self.instance.abs_diff(self, current, self.kmer_array(kmer))

    def to_dict(self, kmer_str=True, params=False):
        if kmer_str:
            d = {"kmer" : self.KMER_STRS}
        else:
            d = {"kmer" : self.KMERS}
        d.update(self._base)
        d.update(self._extra)
        if params:
            d.update(self.params_to_dict("_"))
        return d

    def params_to_dict(self, prefix=""):
        d = dict()
        for name in PARAM_TYPES.keys():
            val = getattr(self.PRMS, name)
            if not isinstance(val, str) or len(val) > 0:
                d[f"{prefix}{name}"] = val 
        return d

    def to_df(self, kmer_str=True):
        return pd.DataFrame({
            key : vals for key,vals in self.to_dict(kmer_str).items()
            if len(vals) > 0
        })

    def to_tsv(self, out=None):
        return self.to_df().to_csv(out, sep="\t", index=False)

    def to_npz(self, fname):
        np.savez_compressed(fname, **self.to_dict(params=True))

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
        d = self.to_dict(params=True)
        return d
    
    def __repr__(self):
        ret = "<PoreModel mean=%.3f stdv=%.3f>\n" % (self.model_mean, self.model_stdv)
        ret += str(self.to_df())
        return ret[:-1]
    
    FILE_LOADERS = {
        "npz" : _vals_from_npz,
        "h5" : _vals_from_hdf5,
        "hdf5" : _vals_from_hdf5,
    }

PoreModel._init_presets()
