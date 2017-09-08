import sys
import numpy as np

REF_START = 3912154
#REF_START = 4641647 

truth_in = open(sys.argv[1])
align_in = open(sys.argv[2])
#align_in = sys.stdin
nofwd = len(sys.argv) > 3 and sys.argv[3] == "nofwd"


event_locs = dict()

truth_in.readline()
for line in truth_in:
    event, ref = map(int, line.split()[1:3])
    event_locs[event] = REF_START - ref + 1

tp_seeds = set()
fp_seeds = set()
max_probs = dict()
max_probs_tp = dict()

tp_probs = list()
fp_probs = list()

unique_ends = set()

not_in = 0

for line in align_in:
    strand, seed_range, ref_range, prob = line.strip().split()
    seed_st, seed_en = map(int, seed_range.split("-"))
    ref_en, ref_st = map(int, ref_range.split("-"))
    prob = float(prob)

    unique_ends.add(seed_en)

    if nofwd and strand == "fwd":
        continue

    if seed_en in event_locs:

        is_max = prob > max_probs.get(seed_en, -1000)

        if is_max:
            max_probs[seed_en] = prob

        if event_locs[seed_en] in range(ref_en-5, ref_en+5):

            if not seed_en in tp_seeds:
                tp_probs.append(prob)

            tp_seeds.add(seed_en)
            
            if is_max:
                max_probs_tp[seed_en] = True

        else:

            if not seed_en in fp_seeds:
                fp_probs.append(prob)

            fp_seeds.add(seed_en)

            if is_max:
                max_probs_tp[seed_en] = False

    else:
        not_in += 1 


#for seed in sorted(max_probs.keys()):
#    print "%s\t%d\t%d\t%d" % (seed, fp_seeds.get(seed, 0), seed in tp_seeds, not max_probs_fp[seed])
total = len(tp_seeds.union(fp_seeds))
unique_tp = len(tp_seeds.difference(fp_seeds))
max_tp = len(filter(lambda t: max_probs_tp[t], max_probs_tp))
print "Total: %d" % (total)
print "TP: %d (%.2f%%)" % (len(tp_seeds), 100.0*len(tp_seeds)/total)
print "FP: %d (%.2f%%)" % (len(fp_seeds), 100.0*len(fp_seeds)/total)
print "max prob TP: %d (%.2f%%)" % (max_tp, 100.0*max_tp/total)
print "unique TP: %d (%.2f%%)" % (unique_tp, 100.0*unique_tp/total)



            
