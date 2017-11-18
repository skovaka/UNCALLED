import h5py
import sys
import numpy as np
import argparse
import os

def get_events(filename):
    events = list()

    hdf = h5py.File(os.path.abspath(filename),'r')

    readnum = filename.split('_')[-2]
    readname='R'+readnum[1:4]+'_'+readnum[4:]

    #events = hdf['Analyses']['EventDetection_000']['Reads'][readname]['Events']
    reads = hdf['Analyses']['EventDetection_000']['Reads']
    
    for r in reads:
        sys.stderr.write("Reading '%s'\n" % r)
        for start, length, mean, stdv in reads[r]['Events']:
            yield (mean, stdv, length)

if __name__ == "__main__":
    i = 0
    for e in get_events(sys.argv[1]):
        print "%d\t%.4f\t%.4f\t%d" % ((i,)+e)
        i += 1
