import sys
import time
import glob
import numpy as np
import pandas as pd
import os
import re
from collections import Sequence

from _uncalled import _PoreModel
from . import nt

class PoreModel(_PoreModel):
    
    def __getitem__(self, kmer):
        return self.means[nt.kmer_array(kmer)]

    #def to_pandas(self):
    #    kmers = np.arange(nt.K)
    #    return pd.DataFrame({"means" : 
