import sys, os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import re
import time
from matplotlib.ticker import NullFormatter, FuncFormatter
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib import widgets
import matplotlib
import scipy.stats
import types
import pandas as pd

from ..sigproc import ProcRead
from ..config import Config, ParamGroup, Opt
from ..index import BWA_OPTS
from ..fast5 import Fast5Reader
from .dtw import Track, ref_coords
from .align import GuidedDTW, BcFast5Aln
from .dotplot import Dotplot
from _uncalled import nt

#BWA_OPTS + 
OPTS = (
    Opt("ref_bounds", "align", type=ref_coords),
    Opt("track_in", type=str),
    Opt(("-m", "--pore-model"), "mapper", default=None),
    Opt("--rna", fn="set_r94_rna"),
    Opt(("-f", "--full-overlap"), "browser", action="store_true"),
    Opt(("-o", "--outfile"), type=str, default=None),
)

def main(conf):
    """Outputs summary statistics comparing aligned current distibutions to the pore model"""

    track = Track(conf.track_in, conf=conf)

    desc = scipy.stats.describe(track[Track.PA_LAYER])#, nan_policy='omit')
    print(desc)
