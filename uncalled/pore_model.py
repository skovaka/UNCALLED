import sys
import time
import glob
import numpy as np
import pandas as pd
import os
import re
from collections import Sequence

from _uncalled import _PoreModel
from . import nt, config

class PoreModel(_PoreModel):

    def __init__(self, name=None, reverse=None, complement=None, df=None, prms=config.DEFAULTS.pore_model):
        if name is not None: prms.name = name
        if reverse is not None: prms.reverse = reverse
        if complement is not None: prms.complement = complement

        if df is not None:
            _PoreModel.__init__(self, np.ravel(df[["mean","stdv"]]), prms.reverse, prms.reverse)
        else:
            _PoreModel.__init__(self, prms)


    def __getitem__(self, kmer):
        return self.means[nt.kmer_array(kmer)]

    def to_df(self):
        return pd.DataFrame({
            "kmer" : nt.KMER_STRS, 
            "mean" : self.means, 
            "stdv" : self.stdvs
        })
