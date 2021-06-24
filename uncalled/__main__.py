"""Utility for Nanopore Read ALignment to Large Expanses of DNA"""

import sys
import numpy as np
import pandas as pd
from . import index, map, realtime, sim, pafstats, dtw, config

SUBCMDS = [
    index, 
    map, 
    realtime, 
    sim, 
    pafstats,
    dtw
]

def main():
    parser = config.ArgParser(SUBCMDS, __doc__)

    cmd, conf = parser.parse_args()

    if cmd is not None:
        ret = cmd(conf=conf)

        if isinstance(ret, pd.DataFrame):
            ret.to_csv(sys.stdout, sep="\t")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
