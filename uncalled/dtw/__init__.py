"""DTW alignment and analysis methods"""

from .dtw import *
from . import align, dotplot, browser, convert #refstats #compare, 

SUBCMDS = [
    align, dotplot, browser, convert, 
    (compare, COMPARE_OPTS), 
    (method_compare, METHOD_COMPARE_OPTS),
    (refstats, REFSTATS_OPTS),
]
