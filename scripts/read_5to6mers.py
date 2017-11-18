import h5py
import sys
import numpy as np
import argparse
import os
from kmer_model import *

def parse_events(model, events):
    events = [(e[4].decode('UTF-8'), e[0], e[2], LEN_FACTOR*e[3]) for e in events]
    
    events_raw = [(mean, stdv, length) for kmer, mean, stdv, length in events]
    
    scale, shift = model.get_norm_params(events_raw)

    #new_means = [(mean - shift) / scale for mean, _, _ in events_raw]
    new_means = [mean for mean, _, _ in events_raw]

    return [(e[0], mean, e[2], e[3]) for mean, e in zip(new_means, events)]

MATCH = 0
STAY = 1
SKIP = 2

LEN_FACTOR = 4000

if __name__ == "__main__":
    i = 0
    model = KmerModel(sys.argv[1])
    filename = sys.argv[2]

    print_skip = len(sys.argv) > 3 and sys.argv[3] == "printskip"


    hdf = h5py.File(filename,'r')

    events = parse_events(model, hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events'])

    event_info = list()


    i = 0
    while i < len(events):

        k1, m1, s1, l1 = events[i] 

        homopoly = True
        for c in range(0, len(k1)-1):
            if k1[c] != k1[c+1]:
                homopoly = False
                break

        #Find next non-stay
        j = i + 1
        while j < len(events) and events[i][0] == events[j][0]:
            j += 1

        if j == len(events):
            break

        #kmer, mean, stdv, length
        k2, m2, s2, l2 = events[j] 

        if k1[1:] == k2[:-1]:

            if homopoly and j - i > 1:
                k3 = k1 + k1[0]
            else:
                k3 = k1 + k2[-1]

            t = 'MATCH' if j - i == 1 else 'STAY'

            for l in range(i, j):
                print ("%s\t%s\t%.4f\t%.4f\t%d" % ((k3, t) + events[l][1:]))

        elif k1[2:] == k2[:-2]:
            k3 = k1 + k2[-2]

            for l in range(i, j):
                print ("%s\t%s\t%.4f\t%.4f\t%d" % ((k3, t) + events[l][1:]))
            skipped = k1[1:] + k2[-2:]

            
            model_stats = model.get_norm_vals(model.kmer_to_i(skipped))

            if print_skip:
                print ("%s\tSKIP\t%.4f\t%.4f\t%.4f\t%.4f" % ((skipped, ) + model_stats))


        i = j




