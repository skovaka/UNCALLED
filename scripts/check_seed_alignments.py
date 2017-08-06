import sys
import numpy as np

REF_START = 3912154

truth_in = open(sys.argv[1])
align_in = open(sys.argv[2])

event_locs = dict()

truth_in.readline()
for line in truth_in:
    event, ref = map(int, line.split()[1:3])
    event_locs[event] = REF_START - ref + 1

tp_seeds = set()
fp_seeds = dict()
max_probs = dict()
max_probs_fp = dict()

tp_probs = list()
fp_probs = list()

unique_ends = set()

for line in align_in:
    strand, seed_range, ref_range, prob = line.strip().split()
    seed_st, seed_en = map(int, seed_range.split("-"))
    ref_st, ref_en = map(int, ref_range.split("-"))
    prob = float(prob)

    unique_ends.add(seed_en)

    #if strand == "fwd":
    #    continue

    if seed_en in event_locs:

        is_max = prob > max_probs.get(seed_en, -1000)
        if is_max:

            max_probs[seed_en] = prob

        if event_locs[seed_en] in range(ref_en-5, ref_en+5):
            if not seed_en in tp_seeds:
                tp_probs.append(prob)
            tp_seeds.add(seed_en)
            
            if is_max:
                max_probs_fp[seed_en] = False

        else:

            if is_max:
                max_probs_fp[seed_en] = True

            if seed_en in fp_seeds:
                fp_seeds[seed_en] += 1
            else:
                fp_probs.append(prob)
                fp_seeds[seed_en] = 1


#for seed in sorted(max_probs.keys()):
#    print "%s\t%d\t%d\t%d" % (seed, fp_seeds.get(seed, 0), seed in tp_seeds, not max_probs_fp[seed])
print len(tp_seeds)
print len(fp_seeds)
print len(unique_ends)
print len(filter(lambda t: max_probs_fp[t], max_probs_fp))
print
print len(filter(lambda p: p <= -3.3, tp_probs)), len(tp_probs)
print len(filter(lambda p: p <= -3.3, fp_probs)), len(fp_probs)



            
