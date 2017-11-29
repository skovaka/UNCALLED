import h5py
import sys
import numpy as np
import argparse
import os

def get_events(filename):
    events = list()

    hdf = h5py.File(os.path.abspath(filename),'r')

    reads = hdf['Analyses']['EventDetection_000']['Reads']
    
    for r in reads:
        sys.stderr.write("Reading '%s'\n" % r)
        for start, length, mean, stdv in reads[r]['Events']:
            yield (length, mean, stdv, start)

if __name__ == "__main__":
    for e in get_events(sys.argv[1]):
        print ("%d\t%.4f\t%.4f\t%d" % e)
