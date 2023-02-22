import pandas as pd

from _uncalled import _AlnCoords, _AlnCoordsDF, _DtwDF, _AlnDF

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

class DtwDF(DataFrameHelper):
    def __init__(self, *args):
        self._init(_DtwDF, *args)

class AlnCoordsDF(_AlnCoordsDF, DataFrameHelper):
    def __init__(self, *args):
        self._init(_AlnCoordsDF, *args)

