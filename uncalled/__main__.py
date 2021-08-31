"""Utility for Nanopore Read ALignment to Large Expanses of DNA"""

import sys
import numpy as np
import pandas as pd
from .argparse import ArgParser

#from . import index, map, realtime, sim, pafstats, dtw
from . import index, pafstats
from .rt import realtime, map, sim
from .dtw import dtw, convert
from .vis import browser, dotplot, sigplot
from .stats import _refstats, dtwstats

SUBCMDS = [
    index, 
    realtime, map, sim, pafstats,
    dtw, convert,
    browser, dotplot, sigplot,
    _refstats#, dtwstats
]
#    (compare, COMPARE_OPTS), 
#    (method_compare, METHOD_COMPARE_OPTS),
#    (refstats, REFSTATS_OPTS),

#commands = [
#    ("General", [
#        ("index", "Builds an UNCALLED index from a FASTA reference"),
#    ],("Real-Time Mapping", [
#        ("realtime", "Perform real-time targeted (ReadUntil) sequencing"),
#        ("map", "Map fast5 files to a DNA reference"),
#        ("sim", "Simulate real-time targeted sequencing"),
#        ("pafstats", "Computes speed and accuracy of UNCALLED mappings"),
#    ]),("DTW Track Generation", [
#        ("dtw", "Aligns reads to a reference guided by basecalled alignments"),
#        ("convert", "Converts tombo or nanopolish alignments to uncalled DTW track format"),
#    ]),("DTW Visualization", [
#        ("browser", "Plot, analyze, and compare DTW alignment tracks interactively"),
#        ("dotplot", "Plot dotplots of DTW alignments"),
#        ("sigplot", "Plot nanopore signal annotated based on DTW alignments"),
#    ]),("DTW Analysis", [
#        ("refstats", "Compute statistics for each reference coordinate based on a DTW track or a pair of tracks")
#    ])
#)]


#SUBCMDS = [
#    dtw, dotplot, sigplot, browser, convert, 
#    (compare, COMPARE_OPTS), 
#    (method_compare, METHOD_COMPARE_OPTS),
#    (refstats, REFSTATS_OPTS),
#]

def main():
    parser = ArgParser(SUBCMDS, __doc__)

    cmd, conf = parser.parse_args()

    if cmd is not None:
        ret = cmd(conf=conf)

        if isinstance(ret, pd.DataFrame):
            ret.to_csv(sys.stdout, sep="\t")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
