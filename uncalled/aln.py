import pandas as pd
import numpy as np

import _uncalled
from _uncalled import _AlnCoords, _AlnCoordsDF, _AlnDF#, _Alignment

from collections import namedtuple
_Layer = namedtuple("_Layer", ["dtype", "label", "default"], defaults=[None,None,False])

ALN_LAYERS = {
    "start" : _Layer("Int32", "Sample Start", True),
    "length" : _Layer("Int32", "Sample Length", True),
    "current" : _Layer(np.float32, "Current (pA)", True),
    "current_sd" : _Layer(np.float32, "Stdv (pA)", True),
    "start_sec" : _Layer(np.float32, "Time (s)", True),
    "length_sec" : _Layer(np.float32, "Dwell (s)", True),
    "length_ms" : _Layer(np.float32, "Dwell (ms)", True),
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

    def __init__(self, seq, name, offset):
        self.instance = seq
        self.name = name
        self.offset = offset

    @property
    def is_flipped(self):
        return self.coords.start < 0

    @property
    def mpos(self):
        return self.coords.expand().to_numpy()

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
        return pd.DataFrame(cols)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"Sequence has no attribute '{name}'")
        return self.instance.__getattribute__(name)

class AlnDF:
    sample_rate = 4000

    def __init__(self, seq, start=None, length=None, current=None, current_sd=None):
        self.seq = seq
        if isinstance(start, _AlnDF):
            self.instance = start
            return

        if current is None:
            current = np.zeros(0, np.float32)
        if current_sd is None:
            current_sd = np.zeros(0, np.float32)

        self.instance = _AlnDF(self.seq.coords, start, length, current, current_sd)

        self._extra = dict()

    @property
    def start(self):
        return self.samples.starts.to_numpy()

    @property
    def end(self):
        return self.samples.ends.to_numpy()

    @property
    def length(self):
        return self.samples.lengths.to_numpy()

    @property
    def start_sec(self):
        return self.start / self.sample_rate

    @property
    def length_sec(self):
        return self.length / self.sample_rate

    @property
    def middle(self):
        return (self.start + self.length / 2).astype(np.float32)

    @property
    def middle_sec(self):
        return self.middle / self.sample_rate

    @property
    def model_diff(self):
        return np.array(self.current) - self.seq.current

    @property
    def model_diff(self):
        return np.array(self.current) - self.seq.current

    @property
    def length_ms(self):
        return 1000* self.length / AlnDF.sample_rate

    def to_pandas(self, layers=None):
        if layers is None:
            df = pd.DataFrame({
                #"mpos" : self.index.expand(),
                "start" : self.samples.starts,
                "length" : self.samples.lengths,
                "current" : self.current,
                "current_sd" : self.current_sd
            })#.set_index("mpos")
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

        return df#.set_index(index)

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

        return df#.set_index(index)

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"AlnDF has no attribute '{name}'")
        return self.instance.__getattribute__(name)

class Alignment:
    def __init__(self, aln_id, read, seq, sam=None):
        self.id = aln_id
        self.seq = seq
        self.sam = sam

        if isinstance(read, str):
            read_id = read
            self.read = None
        else:
            read_id = read.id
            self.read = read

        Super = getattr(_uncalled, f"_AlignmentK{seq.K}", None)
        if Super is None:
            raise ValueError(f"Invalid k-mer length {seq.K}")
        self.instance = Super(read_id, seq.instance)

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

        for name in layers.unique(0):
            _,group_layers = layers.get_loc_level(name)
            group = getattr(self, name, [])
            if len(group) > 0:
                vals[name] = group.to_pandas(group_layers)

        df = pd.concat(vals, axis=1, names=["group", "layer"])

        df["aln","id"] = self.id
        
        if index is not None:
            df.set_index(list(index), inplace=True)
            if join_index:
                df.index.names = [".".join(idx) for idx in index]

        return df#.set_index(index)

    Attrs = namedtuple("Attrs", [
        "id", "read_id", "ref_name", "ref_start", "ref_end", 
        "fwd", "sample_start", "sample_end", "coords"
    ])

    #@property
    def attrs(self):
        return self.Attrs(
            self.id, self.read_id, self.seq.name, self.seq.coords.start, self.seq.coords.end,
            self.seq.fwd, self.moves.samples.start, self.moves.samples.end, self.seq.coords
        )

