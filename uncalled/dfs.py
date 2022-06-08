import pandas as pd

from _uncalled import _AlnCoords

class DataFrameHelper:
    def _init(self, cls, df):
        cls.__init__(self, df[cls.columns].to_records(index=False))

    def to_df(self):
        return pd.DataFrame.from_records(self.to_numpy())

class AlnCoords(_AlnCoords, DataFrameHelper):
    def __init__(self, *args):
        self._init(_AlnCoords, *args)
