import pandas as pd

from _uncalled import _AlnCoords, _AlnCoordsDF

class RecArrayHelper:
    def _init(self, cls, df):
        cls.__init__(self, df[cls.columns].to_records(index=False))

    def to_df(self):
        return pd.DataFrame.from_records(self.to_numpy())

class AlnCoords(_AlnCoords, RecArrayHelper):
    def __init__(self, *args):
        self._init(_AlnCoords, *args)

class DataFrameHelper:
    def _init(self, cls, *args):
        if isinstance(args[0], pd.DataFrame):
            df = args[0]
            args = [df[col].to_numpy(copy=True) for col in self.names]
        cls.__init__(self, *args)

    def __getitem__(self, name):
        return self.__getattribute__(name)

    def to_df(self):
        return pd.DataFrame(
            {name : self[name].to_numpy() for name in self.names},
        )

class AlnCoordsDF(_AlnCoordsDF, DataFrameHelper):
    def __init__(self, *args):
        self._init(_AlnCoordsDF, *args)

