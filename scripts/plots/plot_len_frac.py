import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
import re
from glob import glob

param_names = {'m' : "model",
               'r' : "reference",
               'l' : "read_list",
               'o' : "out_prefix",
               't' : "tally_dist",
               's' : "min_seed_len",
               'a' : "anchor_len",
               'i' : "max_ignores",
               'k' : "max_skips",
               'y' : "stay_frac",
               'E' : "extend_evpr",
               'A' : "anchor_evpr",
               'S' : "seed_pr",
               'Y' : "stay_pr"}



filenames = sys.argv[1:]


for fname in filenames:

    lengths = list()
    fracs = list()

    infile = open(fname)

    while True:
        infile.readline()
        line = infile.readline()
        print ("f", line.strip())

        tabs = line.split()
        if len(tabs) == 0:
            break
        event_ct = int(tabs[-3])

        top_fwd = infile.readline().strip().split()
        print ("s", top_fwd)
        top_rev = None
        next_fwd = None
        next_rev = None

        eof = True

        if len(top_fwd) > 0 and top_fwd[0][0] != "=":
            print("AH")

            for line in infile:
                tabs = line.strip().split()
                eof = False
                if tabs[0] == "fwd":
                    if next_fwd == None and float(tabs[-1]) < 0.9:
                        next_fwd = tabs
                else:
                    if tabs[0] == "rev":
                        top_rev = tabs
                    break

            if eof:
                break
            
            if top_rev != None:
                for line in infile:
                    tabs = line.strip().split()
                    if tabs[0] == "rev":
                        if next_rev == None and float(tabs[-1]) < 0.9:
                            next_rev = tabs
                    else:
                        break
        else:
            tabs = top_fwd
            top_fwd = None

        if top_fwd != None and top_rev != None:
            if int(top_fwd[1]) > int(top_rev[1]):
                best = top_fwd
                next_best = next_fwd
            else:
                best = top_rev
                next_best = next_rev
        elif top_fwd:
            best = top_fwd
            next_best = next_fwd
        elif top_rev:
            best = top_rev
            next_best = next_rev
        else:
            best = None

        print(tabs)

        if best != None and event_ct < 10000:
            aln_length = int(best[1])
            st, en = map(int, best[2].split("-"))
            fracs.append(aln_length / (en - st + 1))
            lengths.append(event_ct)


plt.scatter(lengths, fracs, s=10, c=(0, 0, 1), alpha=0.4, linewidths=0)

plt.show()
