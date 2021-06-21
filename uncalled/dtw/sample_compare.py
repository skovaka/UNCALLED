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
from _uncalled import BwaIndex, nt

#BWA_OPTS + 
OPTS = (
    Opt("ref_bounds", "align", type=ref_coords),
    Opt("track_a", "browser"),
    Opt("track_b", "browser"),
    Opt(("-f", "--full-overlap"), "browser", action="store_true"),
    Opt(("-o", "--outfile"), type=str, default=None),
)

def main(conf):
    """Outputs a TSV file conaining Kolmogorov–Smirnov test statistics comparing the current and dwell time of two alignment tracks"""

    track_a = Track(conf.browser.track_a, conf=conf)
    mat_a = track_a.get_matrix()
    conf.align.mm2_paf = None

    track_b = Track(conf.browser.track_b, conf=conf)
    mat_b = track_b.get_matrix()
    
    # , dtype=[   print(ks.dtype)
    ks = pd.DataFrame(mat_a.calc_ks(mat_b).T, columns=[TrackMatrix.LAYER_META[l][0] for l in TrackMatrix.KS_LAYERS])
    
    csv = ks.to_csv(conf.outfile, sep="\t", index=False)
    if conf.outfile is None:
        print(csv)

