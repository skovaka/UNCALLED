"""DTW alignment and analysis methods"""

from .read_aln import *
from .track import *
from .browser import Browser
from .dotplot import Dotplot
from . import align, convert #refstats #compare, 

SUBCMDS = [
    align, dotplot, browser, convert, 
    (compare, COMPARE_OPTS), 
    (method_compare, METHOD_COMPARE_OPTS),
    (refstats, REFSTATS_OPTS),
]

