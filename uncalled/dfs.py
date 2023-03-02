import pandas as pd
import numpy as np

import _uncalled
from _uncalled import _AlnCoords, _AlnCoordsDF, _AlnDF#, _Alignment

class AlnDF:
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
        return self.samples.starts

    @property
    def end(self):
        return self.samples.ends

    @property
    def length(self):
        return self.samples.lengths

    @property
    def model_diff(self):
        return self.current - self.seq.current

    def to_pandas(self, index="ref"):
        df = pd.DataFrame({
            #"mref" : self.index.expand(),
            "start" : self.samples.starts,
            "length" : self.samples.lengths,
            "current" : self.current,
            "current_sd" : self.current_sd
        })#.set_index("mref")
        
        idx = self.index.expand()
        if index == "ref" and idx[0] < 0:
            idx = -idx-1
        elif index != "mref":
            raise ValueError(f"Unknown index column: {index}")
        df[index] = idx

        return df.set_index(index)

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

    def to_pandas(self, index="ref"):
        df = pd.DataFrame({
            "mean_ref_dist" : self.instance.dist, 
            "jaccard" : self.instance.jaccard,
        })
        
        idx = self.index.expand()
        if index == "ref" and idx[0] < 0:
            idx = -idx-1
        elif index != "mref":
            raise ValueError(f"Unknown index column: {index}")
        df[index] = idx

        return df.set_index(index)

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"AlnDF has no attribute '{name}'")
        return self.instance.__getattribute__(name)

class Alignment:
    def __init__(self, read, seq, sam=None):
        if isinstance(read, str):
            read_id = read
            self.read = None
        else:
            read_id = read.id
            self.read = read

        self.sam = sam

        Super = getattr(_uncalled, f"_AlignmentK{seq.K}", None)
        if Super is None:
            raise ValueError(f"Invalid k-mer length {seq.K}")
        self.instance = Super(read_id, seq)

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

    def to_pandas(self, index="ref", layers=None):
        vals = dict()

        if layers is None:
            for tgt,group in [("dtw","dtw"), ("bcaln","moves"), ("bc_cmp","mvcmp")]:
                g = getattr(self, group, [])
                if len(g) > 0:
                    vals[tgt] = g.to_pandas(index)

        #if len(self.dtw) > 0:
        #    vals["dtw"] = self.dtw.to_pandas(index)

        #if len(self.moves) > 0:
        #    vals["bcaln"] = self.moves.to_pandas(index)

        #if len(self.mvcmp) > 0:
        #    vals["bc_cmp"] = self.mvcmp.to_pandas(index)

        df = pd.concat(vals, axis=1, names=["group", "layer"])

        df["dtw", "kmer"] = self.seq.kmer
        return df
