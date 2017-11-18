import h5py
import sys
import numpy as np
import argparse
import os
from kmer_model import *

MATCH = 0
STAY = 1
SKIP = 2

LEN_FACTOR = 4000

if __name__ == "__main__":
    i = 0
    fast5_dir = sys.argv[1]
    seed_lens = list(map(int, sys.argv[2:]))

    for filename in os.listdir(fast5_dir):
        path = os.path.join(os.path.abspath(fast5_dir), filename)

        hdf = h5py.File(path,'r')

        #mean, stdv, start, length, kmer
        events = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']

        bc_fields = events.dtype.fields
        
        event_info = list()

        p_kmer = None
        p_length = None

        for e in events:
            kmer = e[4].decode('UTF-8')
            length = e[3]*LEN_FACTOR

            if not p_kmer or p_kmer[1:] == kmer[:-1]:
                event_info.append(MATCH)

            elif p_kmer == kmer:
                event_info.append(STAY)

            elif p_kmer[2:] == kmer[:-2]:
                event_info.append(SKIP)

            p_kmer = kmer

        skip_counts = [0]*len(seed_lens)

        for i in range(0, len(event_info)):

            for j in range(0, len(seed_lens)):
                if (i+seed_lens[j] <= len(event_info) and 
                    SKIP in event_info[i:i+seed_lens[j]]):

                    skip_counts[j] += 1


        out = "%d" % len(event_info)
        for sc in skip_counts:
            out += "\t%.1f%%" % (100.0 * sc / len(event_info))

        print (out)

        sys.stdout.flush()

    






