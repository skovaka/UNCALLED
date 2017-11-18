import h5py
import sys
import numpy as np
import argparse
import os

def get_events(filename):
    events = list()

    hdf = h5py.File(os.path.abspath(filename),'r')

    reads = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']
    
    for r in reads:
        sys.stderr.write("Reading '%s'\n" % r)
        for start, length, mean, stdv in reads[r]['Events']:
            yield (length, mean, stdv, start)

if __name__ == "__main__":
    hdf = h5py.File(os.path.abspath(sys.argv[1]),'r')

    events = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']
    
    for e in events:
        length = e['length']
        if int(length) == 0:
            length = int(4000*length)
        print "%s\t%.2f\t%.2f\t%d" % (e['model_state'], e['mean'], e['stdv'], length)
