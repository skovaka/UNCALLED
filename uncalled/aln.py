import pandas as pd
import numpy as np

import _uncalled

from collections import namedtuple
_Layer = namedtuple("_Layer", ["dtype", "label", "default"], defaults=[None,None,False])

ALN_LAYERS = {
    "start" : _Layer("Int32", "Sample Start", True),
    "length" : _Layer("Int32", "Sample Length", True),
    "current" : _Layer(np.float32, "Current (pA)", True),
    "current_sd" : _Layer(np.float32, "Stdv (pA)", True),
    "start_sec" : _Layer(np.float32, "Time (s)", True),
    "length_sec" : _Layer(np.float32, "Dwell (s)", True),
    "dwell" : _Layer(np.float32, "Dwell (ms)", True),
    "end" : _Layer("Int32", "Sample End", False),
    "middle" : _Layer(np.float32, "Sample Middle", False),
    "middle_sec" : _Layer(np.float32, "Sample Middle", False),
    "model_diff" : _Layer(np.float32, "Model pA Diff.", True),
    "abs_diff" : _Layer(np.float32, "Abs. Model Diff.", False),
}

LAYERS = {
    "seq" : {
        "mpos" : _Layer("Int64", "Mirror Ref.", True),
        "pos" : _Layer("Int64", "Reference Coord.", False),
        "pac" : _Layer("Int64", "Packed Ref. Coord.", False),
        "kmer" : _Layer("Int32", "Kmer", True),
        "current" : _Layer(str, "Model Mean (pA)", False),
        "name" : _Layer(str, "Reference Name", False),
        "fwd" : _Layer(bool, "Forward", False),
        "strand" : _Layer(str, "Strand", False),
        "base" : _Layer(str, "Base", False),
    }, "aln" : {
        "id" : _Layer("Int64", "Aln. ID"),
    }, "dtw" : ALN_LAYERS, "moves" : ALN_LAYERS,
    #}, "moves" : {
    #    "start" : _Layer("Int32", "BC Sample Start", True),
    #    "length" : _Layer("Int32", "BC Sample Length", True),
    #    "end" : _Layer("Int32", "Sample End", False),
    #    "middle" : _Layer(np.float32, "Sample Middle", False),
    #    "indel" : _Layer("Int32", "Basecalled Alignment Indel", False),
    "cmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment V", True),
        "aln_b" : _Layer("Int32", "Compare alignment A", True),
        "group_b" : _Layer(str, "Compare type", False),
        "jaccard" : _Layer(np.float32, "Jaccard Distance", True),
        "dist" : _Layer(np.float32, "Mean Ref. Distance", True),
    }, "mvcmp" : {
        "jaccard" : _Layer(np.float32, "Jaccard Distance", True),
        "dist" : _Layer(np.float32, "Mean Ref. Distance", True),
    }
}

LAYER_META = pd.concat([
    pd.concat({
        group : pd.DataFrame(layers, index=_Layer._fields).transpose()
    }, names=("group","layer"))  
    for group, layers in LAYERS.items()
])

DEFAULT_LAYERS = LAYER_META.index[LAYER_META["default"]]

def parse_layer(layer):
    if isinstance(layer, str):
        spl = layer.split(".")
    elif isinstance(layer, tuple):
        spl = layer
    else:
        raise ValueError("Layer must be string or tuple")

    if len(spl) == 2:
        group,layer = spl
    elif len(spl) == 1:
        if layer in LAYER_META.index.get_level_values(0):
            group = layer
            _,layers = DEFAULT_LAYERS.get_loc_level(layer)
            for layer in layers:
                yield (group, layer)
            return
        group = "dtw"
    else:
        raise ValueError("Invalid layer: \"{layer}\"")

    if not (group, layer) in LAYER_META.index:
        raise ValueError(f"Invalid layer \"{group}.{layer}\"")

    yield (group, layer)

def parse_layers(layers):
    if layers is None:
        return pd.Index([])

    if isinstance(layers, str):
        layers = layers.split(",")

    ret = list()

    for layerstr in layers:
        for layer in parse_layer(layerstr):
            ret.append(layer)

    return pd.Index(ret).unique()

class Sequence:
    LAYERS = {"pos", "mpos", "pac", "name", "fwd", "strand", "kmer", "current"}
    CONST_LAYERS = {"name", "fwd", "strand"}
    DEFAULT_LAYERS = ["pos", "kmer"]

    def __init__(self, seq, offset):
        self.instance = seq
        self.offset = offset
        self.index = self.instance.mpos

    @property
    def name(self):
        return self.coord.name

    @property
    def is_flipped(self):
        return self.index.start < 0

    @property
    def mpos(self):
        return self.index.expand().to_numpy()

    @property
    def pos(self):
        if self.is_flipped:
            return -self.mpos-1
        return self.mpos

    @property
    def pac(self):
        return self.offset + self.pos

    @property
    def strand(self):
        return "+" if self.fwd else "-"

    @property
    def fwd(self):
        return self.is_fwd

    def __len__(self):
        return len(self.instance)

    def _iter_layers(self, names):
        ret = list()

    def to_pandas(self, layers=None):
        if layers is None:
            layers = ["seq.pos", "kmer"]

        cols = dict()
        for name in layers:
            val = getattr(self, name)
            if name in self.CONST_LAYERS:
                val = np.full(len(self), val)
            cols[name] = val
        cols["index"] = self.mpos
        return pd.DataFrame(cols).set_index("index")

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"Sequence has no attribute '{name}'")
        return self.instance.__getattribute__(name)

PANDAS_DTYPES = {
    "int16" : pd.Int16Dtype(),
    "int32" : pd.Int32Dtype(),
    "int64" : pd.Int64Dtype(),
    "uint16" : pd.UInt16Dtype(),
    "uint32" : pd.UInt32Dtype(),
    "uint64" : pd.UInt64Dtype(),
}

class AlnDF:
    #sample_rate = 4000

    def __init__(self, seq, start=None, length=None, current=None, current_sd=None):
        self.seq = seq
        if isinstance(start, _uncalled._AlnDF):
            self.instance = start
        else:
            if current is None:
                current = np.zeros(0, np.float32)
            if current_sd is None:
                current_sd = np.zeros(0, np.float32)

            self.instance = _uncalled._AlnDF(self.seq.index, start, length, current, current_sd)

        self.instance.mask(self.na_mask)
        self._extra = dict()

        
    @property
    def na_mask(self):
        return self.samples.mask.to_numpy()

    @property
    def start(self):
        return (self.samples.starts.to_numpy())

    @property
    def end(self):
        return (self.samples.ends.to_numpy())

    @property
    def length(self):
        return (self.samples.lengths.to_numpy())

    @property
    def start_sec(self):
        return (self.start / self.seq.model.PRMS.sample_rate)

    @property
    def length_sec(self):
        return (self.length / self.seq.model.PRMS.sample_rate)

    @property
    def middle(self):
        return (self.start + self.length / 2).astype(np.float32)

    @property
    def middle_sec(self):
        return self.middle / self.seq.model.PRMS.sample_rate

    @property
    def model_diff(self):
        return np.array(self.current) - self.seq.current

    @property
    def model_diff(self):
        return np.array(self.current) - self.seq.current

    @property
    def dwell(self):
        return 1000 * self.length / self.seq.model.PRMS.sample_rate

    def _get_series(self, name):
        vals = getattr(self, name, None)
        if vals is None or len(vals) == 0:
            return None
        ret = pd.Series(vals, copy=True)
        dtype = PANDAS_DTYPES.get(ret.dtype.name, None)
        na = np.nan
        if dtype is not None:
            ret = ret.astype(dtype)
            na = pd.NA
        if len(self.na_mask) > 0:
            ret[~self.na_mask] = na
        return ret

    def to_pandas(self, layers=None):
        if layers is None:
            layers = ["start", "length", "current", "current_sd"]

        df = pd.DataFrame({
            name : self._get_series(name) 
            for name in layers if hasattr(self,name)
        })
        df["index"] = self.index.expand()
        
        #idx = self.index.expand().to_numpy()
        #if index == "ref" and idx[0] < 0:
        #    idx = -idx-1
        #elif index != "mpos":
        #    raise ValueError(f"Unknown index column: {index}")
        #df[index] = idx

        return df.set_index("index")

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"AlnDF has no attribute '{name}'")
        return self.instance.__getattribute__(name)

class CmpDF:
    def __init__(self, seq, instance):
        self.seq = seq
        self.instance = instance

    def to_pandas(self, layers=None):
        if layers is None:
            df = pd.DataFrame({
                "dist" : self.instance.dist, 
                "jaccard" : self.instance.jaccard,
            })
        else:
            df = pd.DataFrame({
                l : getattr(self, l) for l in layers if hasattr(self,l)
            })
        
        #idx = self.index.expand().to_numpy()
        #if index == "ref" and idx[0] < 0:
        #    idx = -idx-1
        #elif index != "mpos":
        #    raise ValueError(f"Unknown index column: {index}")
        #df[index] = idx
        df["index"] = self.index.expand()

        return df.set_index("index")

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"AlnDF has no attribute '{name}'")
        return self.instance.__getattribute__(name)

class Alignment:
    def __init__(self, aln_id, read, seq, sam=None):
        #self.id = aln_id
        self.seq = seq
        self.sam = sam

        if isinstance(read, str):
            read_id = read
            self.read = None
        else:
            read_id = read.id
            self.read = read

        if isinstance(self.seq.model, _uncalled.PoreModelU16):
            Super = _uncalled._AlignmentU16
        elif isinstance(self.seq.model, _uncalled.PoreModelU32):
            Super = _uncalled._AlignmentU32
        else:
            raise ValueError(f"Unknown PoreModel type: {model.instance}")

        #Super = getattr(_uncalled, f"_AlignmentK{seq.K}", None)
        #if Super is None:
        #    raise ValueError(f"Invalid k-mer length {seq.K}")
        self.instance = Super(aln_id, read_id, seq.instance)

        self.dtw = AlnDF(seq, self.instance._dtw)
        self.moves = AlnDF(seq, self.instance._moves)
        self.mvcmp = CmpDF(seq, self.instance._mvcmp)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"Alignment has no attribute '{name}'")
        return self.instance.__getattribute__(name)
        
    def set_dtw(self, df):
        if isinstance(df, AlnDF):
            df = df.instance
        self.instance.set_dtw(df)
        
    def set_moves(self, df):
        if isinstance(df, AlnDF):
            df = df.instance
        self.instance.set_moves(df)

    def to_pandas(self, layers=None, index=None, join_index=True):
        vals = dict()

        layers = parse_layers(layers)
        if index is not None:
            index = parse_layers(index)
            layers = layers.union(index)

        idx = self.seq.mpos#.expand().to_numpy()

        for name in layers.unique(0):
            _,group_layers = layers.get_loc_level(name)
            group = getattr(self, name, [])
            if len(group) > 0:
                vals[name] = group.to_pandas(group_layers).reindex(idx)#.reset_index(drop=True)
                #if idx is None:
                #    idx = vals[name].index
                #else:
                #    idx = vals[name].index.intersection(idx)

        df = pd.concat(vals, axis=1, names=["group", "layer"]).reset_index(drop=True)

        df["aln","id"] = self.id
        
        if index is not None:
            df.set_index(list(index), inplace=True)
            if join_index:
                df.index.names = [".".join(idx) for idx in index]

        return df#.set_index(index)

    Attrs = namedtuple("Attrs", [
        "id", "read_id", "ref_name", "ref_start", "ref_end", 
        "fwd", "sample_start", "sample_end", "coord"
    ])

    #@property
    def attrs(self):
        samp_start = 1000000000
        samp_end = 0
        for df in [self.dtw, self.moves]:
            if not df.empty():
                samp_start = min(samp_start, df.samples.start)
                samp_end = max(samp_end, df.samples.end)

        return self.Attrs(
            self.id, self.read_id, self.seq.coord.name, self.seq.coord.start, self.seq.coord.end,
            self.seq.fwd, samp_start, samp_end, self.seq.coord
        )

