import pandas as pd
import numpy as np

import _uncalled
from _uncalled import _AlnCoords, _AlnCoordsDF, _AlnDF#, _Alignment

class RecArrayHelper:
    def _init(self, cls, df):
        cls.__init__(self, df[cls.columns].to_records(index=False))

    def to_df(self):
        return pd.DataFrame.from_records(self.to_numpy())

class AlnCoords(_AlnCoords, RecArrayHelper):
    def __init__(self, *args):
        self._init(_AlnCoords, *args)

class DataFrameHelper:
    #def _init(self, cls, *args):
    def _init(self, Cls, *args):
        if len(args) == 1 and isinstance(args[0], pd.DataFrame) and Cls is not None:
            df = args[0]
            args = [df[col].to_numpy(copy=True) for col in self.names]
            self.instance = Cls(*args)
        elif len(args) > 1:
            self.instance = Cls(*args)
        elif len(args) == 1:
            self.instance = args[0]

        else:
            raise ValueError(f"Invalid {Cls} dataframe")

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def __getitem__(self, name):
        return getattr(self.instance, name)
        #return self.instance.__getattribute__(name)

    def to_df(self):
        return pd.DataFrame(
            {name : self[name].to_numpy() for name in self.names if len(self[name]) > 0},
        )

class AlnCoordsDF(_AlnCoordsDF, DataFrameHelper):
    def __init__(self, *args):
        self._init(_AlnCoordsDF, *args)

class AlnDF:
    def __init__(self, coords, start=None, length=None, current=None, current_sd=None):
        if isinstance(coords, _AlnDF):
            self.instance = coords
            return

        if current is None:
            current = np.zeros(0, np.float32)
        if current_sd is None:
            current_sd = np.zeros(0, np.float32)
        self.instance = _AlnDF(coords, start, length, current, current_sd)

    def to_pandas(self):
        return pd.DataFrame({
            "mref" : self.index.expand(),
            "start" : self.samples.starts,
            "length" : self.samples.lengths,
            "current" : self.current,
            "current_sd" : self.current_sd
        }).set_index("mref")

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

        self.dtw = AlnDF(self.instance._dtw)
        self.moves = AlnDF(self.instance._moves)

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

    def to_pandas(self, index="ref"):
        vals = {
            "dtw" : self.dtw.to_pandas(), 
            "bcaln" : self.moves.to_pandas(),
        }

        if not self._mvcmp.empty():
            vals["bc_cmp"] = pd.DataFrame({
                    "mref" : self._mvcmp.index.expand(),
                    "mean_ref_dist" : self._mvcmp.dist, "jaccard" : self._mvcmp.jaccard,
            }).set_index(["mref"])

        df = pd.concat(vals, axis=1, names=["group", "layer"])

        df["dtw", "kmer"] = self.seq.kmer
        return df
