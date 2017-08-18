import sys

MATCH_LEN = 16

seeds_in = open(sys.argv[1])

seed_counts = dict()
total_count = 0

for line in seeds_in:
    _, event_range, ref_range, _ = line.strip().split()

    e_st, e_en = map(int, event_range.split("-"))
    seed_len = e_en - e_st + 1
    stay_len = seed_len - MATCH_LEN
    s = (seed_len, stay_len)

    if s in seed_counts:
        seed_counts[s] += 1
    else:
        seed_counts[s] = 1
    total_count += 1

prc_count = 0

for seed_len, stay_len in sorted(seed_counts):
    sc = seed_counts[(seed_len, stay_len)]
    prc_count += sc


    print "%d\t%d\t%.2f\t%.2f\t%d" % (seed_len, stay_len,
                                float(stay_len) / seed_len,
                                float(prc_count) / total_count,
                                sc)


