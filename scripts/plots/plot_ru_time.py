import sys
import matplotlib.pyplot as plt
from scipy.cluster.vq import *
import numpy as np
import os, re

REF_LEN = 4641647 

read_id_re = re.compile("_(ch\d+_read\d+)_")

read_locs_in = open(sys.argv[1])

read_locs = dict()
for line in read_locs_in:
    read_id, loc, length = line.split()
    l = tuple(map(int, (loc, length)))

    if read_id in read_locs:
        read_locs[read_id].append(l)
    else:
        read_locs[read_id] = [l]

good_lens = list()
good_times = list()

bad_lens = list()
bad_times = list()

no_seed_count = 0
no_conf_count = 0
uncluded_count = 0

for line in sys.stdin:
    tabs = line.split()

    if tabs[1] == "FAILED":
        if int(tabs[2]) > 0:
            no_conf_count += 1
        else:
            no_seed_count += 1
        continue

    filename, strand, time, ref_loc, evt_loc, length, ratio = tabs
    print(tabs)
    read_id = read_id_re.search(filename).group(1)
    time = float(time)
    ref_st, ref_en = map(int, ref_loc.split("-"))
    if evt_loc[0] == "-":
        evt_st = 0
        evt_en = int(evt_loc.split("-")[-1])
    else:
        evt_st, evt_en = map(int, evt_loc.split("-"))
    length = int(length)
    ratio = float(ratio)

    if strand == "rev":
        tmp = ref_st
        ref_st = REF_LEN - ref_en
        ref_en = REF_LEN - tmp

    correct = False
    if read_id in read_locs:
        for r_st, r_len in read_locs[read_id]:
            if ref_en > r_st and ref_st < r_st+r_len:
                correct = True
                break



    if correct:
        good_lens.append(4000-evt_st+1)
        good_times.append(time/1000.0)
    else:
        bad_lens.append(4000-evt_st+1)
        bad_times.append(time/1000.0)

    if time > 40000:
        uncluded_count += 1

total_count = float(len(good_times) + len(bad_times) + no_seed_count + no_conf_count)

print (uncluded_count)
print ("Correctly Aligned:\t%d (%.2f%%)" % (len(good_times), 100 * (len(good_times)) / total_count))
print ("Incorrectly Aligned:\t%d (%.2f%%)" % (len(bad_times), 100 * len(bad_times) / total_count))
print ("Not Aligned:\t%d (%.2f%%)" % (no_seed_count + no_conf_count, 100 * (no_seed_count + no_conf_count) / total_count))

plt.scatter(good_lens, good_times, s=40, c='blue', alpha=0.4, linewidths=0)
plt.scatter(bad_lens, bad_times, s=40, c='red', alpha=1, linewidths=0)
plt.axis([0, 4000, 0, 40])
plt.xlabel("Number of Events Processed", fontsize=30)
plt.ylabel("Time (sec)", fontsize=30)
plt.title('Time/Event Count per Read Required to Confidently Align at Least 100bp', fontsize=30)
plt.show()
