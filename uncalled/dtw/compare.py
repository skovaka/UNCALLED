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
