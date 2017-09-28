import h5py
import sys
import numpy as np
import argparse
import os

if __name__ == "__main__":
    i = 0
    filename = sys.argv[1]

    hdf = h5py.File(os.path.abspath(filename),'r')

    readnum = filename.split('_')[-2]
    readname='R'+readnum[1:4]+'_'+readnum[4:]

    #events = hdf['Analyses']['EventDetection_000']['Reads'][readname]['Events']
    reads = hdf['Analyses']['EventDetection_000']['Reads']
    events = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']
    bc_fields = list(events.dtype.fields.keys())
    #print(" ".join(bc_fields))
    p_kmer = None
    for e in events:
        kmer = e[4].decode('UTF-8')
        if not p_kmer or p_kmer[1:] == kmer[:-1]:
            print("match")
        elif p_kmer == kmer:
            print("stay")
        else:
            print("skip")

        p_kmer = kmer
    
