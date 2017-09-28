import sys
from bisect import bisect_left

SUPPORT_THRESH = 8

rev_locs = list() #[ref_prev, read_en, read_prev, sup]
fwd_locs = list()

l = 0

event_count = None
last_out = 0.0

for line in sys.stdin:
    if line[0] == "=":

        if "events" in line:
            event_count = int(line.split()[1])

        continue

    strand, evt_loc, ref_loc, prob = line.strip().split()
    evt_st, evt_en = map(int, evt_loc.split("-"))
    ref_st, ref_en = map(int, ref_loc.split("-"))
    prob = float(prob)
    
    if event_count:
        frac_done = float(event_count - evt_en) / event_count
        if (frac_done - last_out) >= 0.001:
            sys.stderr.write("~ %.1f%% complete\n" % (100*frac_done))
            last_out = frac_done

    if strand == "fwd":
        locs = fwd_locs
    else:
        locs = rev_locs

    new_seed = (ref_en, evt_en, ref_en, 1)

    to_place = -1 #to place ends up >= i
    #print l

    for i in range(0, len(locs)):
        ref_prev, evt_prev, read_en, sup = locs[i]

        higher_sup = to_place < 0 or locs[to_place][3] < sup

        in_range = ref_prev >= ref_en and ref_prev - evt_prev <= ref_en - evt_en

        if higher_sup and in_range:
            to_place = i

        i += 1

    if to_place >= 0 and evt_en != locs[to_place][1]:

        ref_st, sup = locs[to_place][2:]
        new_seed = (ref_en, evt_en, ref_st, sup + 1)

        locs[to_place] = new_seed
    else:
        locs.append(new_seed)

    #else:

    #    locs.append( list(seed_info) )
        #print l, locs[-1]

    l += 1

for ref_en, evt_en, ref_st, sup in sorted(rev_locs, key=lambda l: -l[3]):
    if sup >= SUPPORT_THRESH:
        print "rev", sup, ref_en, ref_st, evt_en

for ref_en, evt_en, ref_st, sup in sorted(fwd_locs, key=lambda l: -l[3]):
    if sup >= SUPPORT_THRESH:
        print "fwd", sup, ref_en, ref_st, evt_en




