import h5py
import sys
import numpy as np
import argparse
import os

def get_len(filename):
    events = list()

    hdf = h5py.File(os.path.abspath(filename),'r')

    reads = hdf['Analyses']['EventDetection_000']['Reads']
    
    if len(reads.keys()) == 1:
        return len(reads[reads.keys()[0]]['Events'])

    lens = list()

    for r in reads:
        lens.append(len(reads[r]['Events']))
    return lens

if __name__ == "__main__":
    d = sys.argv[1]
    for fname in os.listdir(d):
        path = os.path.join(d, fname)
        print "%s\t%s" % (path, get_len(path))
