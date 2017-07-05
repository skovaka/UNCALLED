import sys
from bisect import bisect_left

events_in = open(sys.argv[1])
align_in = open(sys.argv[2])

lvl_ids = dict()
events = list()
for line in events_in:
    tabs = line.strip().split()
    i, mean, stdv = int(tabs[0]), float(tabs[1]), float(tabs[2])
    events.append((mean, stdv, i))

events_rev = list(reversed(events))

prev_evt = -1

align_in.readline()
for line in align_in:
    tabs = line.strip().split()
    ref_contig, ref_pos = tabs[0], int(tabs[1])
    event_mean, event_stdv = float(tabs[6]), float(tabs[7])
    event_kmer = tabs[9]

    rev = tabs[4] in ['c', 't']

    evt = prev_evt+1
    for mean, stdv, j in (events_rev if rev else events)[prev_evt+1:]:
        if abs(event_mean - mean) < 0.01 and abs(event_stdv - stdv) < 0.001:
            break
        evt += 1

    prev_evt = evt

    if rev:
        evt = len(events) - evt

    print ("\t".join(map(str, [evt, ref_pos, event_mean, event_stdv, event_kmer])))
