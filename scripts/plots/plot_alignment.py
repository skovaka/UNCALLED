import numpy as np
import matplotlib.pyplot as plt
import sys

evt_locs = list()
ref_locs = list()
colors = list()

prev_loc = None

st = 0

for line in sys.stdin:
    tabs = line.strip().split()
    evt_loc, ref_loc, mean, stdv, kmer = int(tabs[0]), int(tabs[1]), float(tabs[2]), float(tabs[3]), tabs[4]
    evt_locs.append(evt_loc)
    ref_locs.append(ref_loc)

    if not prev_loc or abs(ref_loc - prev_loc) == 1:
        colors.append("green")
    elif abs(ref_loc - prev_loc) == 0:
        colors.append("blue")
    else:
        colors.append("red")

    prev_loc = ref_loc

plt.scatter(evt_locs, ref_locs, c=colors, s=2)

plt.show()
