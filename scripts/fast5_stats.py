import h5py
import sys
import numpy as np
import argparse
import os
from kmer_model import *

MATCH = 0
STAY = 1
SKIP = 2
IGNORE = 3

def parse_basecalls(model, basecalls):
    basecalls = [(e[4].decode('UTF-8'), e[0], e[2], LEN_FACTOR*e[3]) for e in basecalls]
    
    basecalls_raw = [(mean, stdv, length) for kmer, mean, stdv, length in basecalls]
    
    scale, shift = model.get_norm_params(basecalls_raw)

    #new_means = [(mean - shift) / scale for mean, _, _ in basecalls_raw]
    new_means = [mean for mean, _, _ in basecalls_raw]

    basecalls = [(e[0], mean, e[2], e[3]) for mean, e in zip(new_means, basecalls)]
    
    labeled = list()
        
    i = 0
    while i < len(basecalls):

        k1, m1, s1, l1 = basecalls[i] 

        homopoly = True
        for c in range(0, len(k1)-1):
            if k1[c] != k1[c+1]:
                homopoly = False
                break

        #Find next non-stay
        j = i + 1
        while j < len(basecalls) and basecalls[i][0] == basecalls[j][0]:
            j += 1

        if j == len(basecalls):
            break

        #kmer, mean, stdv, length
        k2, m2, s2, l2 = basecalls[j] 

        if k1[1:] == k2[:-1]:

            if homopoly and j - i > 1:
                k3 = k1 + k1[0]
            else:
                k3 = k1 + k2[-1]

            t = MATCH if j - i == 1 else STAY

            labeled.append((k3, MATCH) + basecalls[i][1:])
            for l in range(i+1, j):
                labeled.append((k3, STAY) + basecalls[l][1:])

        elif k1[2:] == k2[:-2]:
            k3 = k1 + k2[-2]

            labeled.append((k3, MATCH) + basecalls[i][1:])
            for l in range(i+1, j):
                labeled.append((k3, STAY) + basecalls[l][1:])

            skipped = k1[1:] + k2[-2:]

            model_stats = model.get_norm_vals(model.kmer_to_i(skipped))

            labeled.append((k3, SKIP) + basecalls[l][1:])

        i = j

    return labeled

def parse_events(reads):
    events = list()

    for r in reads:
        sys.stderr.write("Reading '%s'\n" % r)
        for start, length, mean, stdv in reads[r]['Events']:
            events.append( (mean, stdv, length) )

    return events

def events_match(e1, e2):
    THRESH = 0.01
    matches = 0
    for i, j in zip(e1, e2):
        if abs(i - j) <= THRESH:
            matches += 1
    return matches >= len(e1) - 1

def print_stats(model, seed_lens, filename):
    hdf = h5py.File(filename,'r')

    basecalls = parse_basecalls(model, hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events'])
    events = parse_events(hdf['Analyses']['EventDetection_000']['Reads'])

    print "== %s ==" % filename

    e = 0

    combined = list()

    for b in range(0, len(basecalls)):

        if basecalls[b][1] == SKIP:
            kmer, t, mean, stdv, length = basecalls[b]
            combined.append( (-1, mean, stdv, length, kmer, t) )
        else:
            while e < len(events) and not events_match(events[e], basecalls[b][2:]):
                mean, stdv, length = events[e]
                combined.append ((e, mean, stdv, length, "NNNNNN", IGNORE))
                e += 1

            k, t = basecalls[b][:2]
            mean, stdv, length = events[e]

            combined.append( (e, mean, stdv, length, k, t) )
            e += 1


    labels = [e[-1] for e in combined]

    for i in range(0, len(combined) - min(seed_lens) + 1):
        for l in seed_lens:
            if i + l < len(combined):
                stays = labels[i:i+l].count(STAY) / float(l)
                skips = labels[i:i+l].count(SKIP) / float(l)
                ignores = labels[i:i+l].count(IGNORE) / float(l)
                #sys.stdout.write("%3.2f %3.2f %3.2f\t" % (stays, skips, ignores))
                #sys.stdout.write("%.5f " % stays)
                sys.stdout.write("%.5f " % ignores)
            else:
                #sys.stdout.write("NA   NA   NA\t")
                sys.stdout.write("NA     ")
        sys.stdout.write("\n")



LEN_FACTOR = 4000

if __name__ == "__main__":
    i = 0
    model = KmerModel(sys.argv[1])

    seed_lens = [8, 16, 24, 32, 48, 64]
    readlist = open(sys.argv[2])

    print "\t".join(map(str, seed_lens))

    for line in readlist:
        print_stats(model, seed_lens, line.strip())



    




