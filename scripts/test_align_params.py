import sys
import numpy as np

EVT_KMER = 0
EVT_ID = 1
REF_LOC = 2
EVT_MEAN = 3
EVT_STDV = 4
LV_MEAN = 5
LV_STDV = 6
SD_MEAN = 7
NM_PROB = 8
IG_PROB = 9
LOG_PROB = 10

MATCH_LEN = int(sys.argv[1])

sys.stdin.readline() #read past header

events = list()

for line in sys.stdin:
    evt = line.strip().split()
    events.append(evt[:1] + list(map(int, evt[1:3])) + list(map(float, evt[3:])))

seed_min = list()
seed_avg = list()

pass_thresh = 0

for i in range(0, len(events) - MATCH_LEN + 1):
    min_log = 0
    avg_log = 0
    stays = 0
    seed_skip = False

    prev_evt_id = None
    prev_ref_id = None

    j = i
    while j - i - stays < MATCH_LEN and j < len(events):
        
        if prev_ref_id != None:
            step = np.abs(events[j][REF_LOC] - prev_ref_id)
        else:
            step = 1

        if step == 0:
            stays += 1
        elif step > 1:
            seed_skip = True
            break

        if prev_evt_id != None and np.abs(events[j][EVT_ID] - prev_evt_id) > 1:
            stays += 1

        min_log = min(min_log, events[j][LOG_PROB])
        avg_log += events[j][LOG_PROB]

        prev_evt_id = events[j][EVT_ID]
        prev_ref_id = events[j][REF_LOC]

        j += 1

    avg_log /= stays + MATCH_LEN

    if not seed_skip and avg_log > -3.75:
        pass_thresh+=1

    if not seed_skip:
        seed_min.append(min_log)
        seed_avg.append(avg_log)
        #print "%d\t%.3f\t%.3f\t%d" % (events[i][EVT_ID], min_log, avg_log, stays)

    #print "%d\t%d\t%.3f\t%.3f\t%d" % (events[i][EVT_ID], not seed_skip, min_log, avg_log, stays)

print len(seed_min)
print min(seed_min)
print min(seed_avg)
print pass_thresh


