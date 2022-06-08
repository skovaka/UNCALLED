import pandas as pd

#from _uncalled import _AlnCoords

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

#class AlnCoords(_AlnCoords, DataFrameHelper):
#    def __init__(self, *args):
#        self._init(_AlnCoords, *args)

