import sys
import numpy as np
import pandas as pd
import textwrap
from .argparse import ArgParser

#from . import index, map, realtime, sim, pafstats, dtw
from . import index, pafstats
from .rt import realtime, map, sim
from .dtw import dtw, convert
from .vis import browser, dotplot, sigplot, trackplot
from .stats import _refstats, _dtwstats, _readstats

SUBCMDS = [
    index, 
    realtime, map, sim, pafstats,
    dtw, convert,
    browser, dotplot, sigplot, trackplot,
    _refstats, _readstats, _dtwstats
]

_help_lines = [
    "Utility for Nanopore Read ALignment to Large Expanses of DNA", "",
    "subcommand options:",
    "General:",
    "\tindex      " + index.main.__doc__.split(".")[0],"",
    "Real-Time Enrichment (ReadUntil/Adaptive Sampling):",
    "\trealtime   " + realtime.main.__doc__,
    "\tmap        " + map.main.__doc__,
    "\tsim        " + sim.main.__doc__,
    "\tpafstats   " + pafstats.main.__doc__, "",
    "Dynamic Time Warping (DTW) Alignment:",
    "\tdtw        " + dtw.main.__doc__,
    "\tconvert    " + convert.__doc__, "",
    "DTW Analysis:",
    "\tdtwstats   " + _dtwstats.main.__doc__,
    "\trefstats   " + _refstats.main.__doc__,
    "\treadstats  " + _readstats.main.__doc__,"",
    "DTW Visualization:",
    "\tbrowser    " + browser.main.__doc__,
    "\tsigplot    " + sigplot.main.__doc__,
    "\tdotplot    " + dotplot.main.__doc__,
]

HELP = "\n".join([
    textwrap.fill(line,
        width=75,
        drop_whitespace=False, 
        replace_whitespace=False, 
        tabsize=2, 
        subsequent_indent=" "*13)
    for line in _help_lines])

def main():
    parser = ArgParser(SUBCMDS, HELP)

    cmd, conf = parser.parse_args()

    if cmd is not None:
        ret = cmd(conf=conf)

        if isinstance(ret, pd.DataFrame):
            ret.to_csv(sys.stdout, sep="\t")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
