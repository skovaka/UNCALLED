"""Rapidly maps raw nanopore signal to dna references"""

import numpy as np
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
        cmd(conf)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
