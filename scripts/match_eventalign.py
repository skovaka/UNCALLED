import sys
from bisect import bisect_left
from extract_events import get_events
from kmer_model import *

fast5_fname = sys.argv[1]
eventalign_in = open(sys.argv[2])

events = list(get_events(fast5_fname))
events_rev = list(reversed(events))

model = KmerModel()

prev_evt = -1

norm_params = model.get_norm_params(events)
print " ".join(["evt_kmer", "evt_id", "ref_loc", "evt_mean", "evt_stdv", "lv_mean", "lv_stdv", "sd_mean", "nm_prob", "ig_prob", "log_prob"])

eventalign_in.readline()

for line in eventalign_in:
    tabs = line.strip().split()
    ref_contig, ref_pos = tabs[0], int(tabs[1])
    evt_mean, evt_stdv = float(tabs[6]), float(tabs[7])
    evt_kmer = tabs[9]

    rev = tabs[4] in ['c', 't']

    evt = prev_evt+1
    for mean, stdv, j in (events_rev if rev else events)[prev_evt+1:]:
        if abs(evt_mean - mean) < 0.01 and abs(evt_stdv - stdv) < 0.001:
            break
        evt += 1

    prev_evt = evt

    if rev:
        evt = len(events) - evt
    
    k_id = model.kmer_to_i(evt_kmer)

    if not k_id:
        lv_mean = lv_stdv = sd_mean = nm_prob = ig_prob = log_prob = 0
    else:
        lv_mean, lv_stdv, sd_mean, ig_lambda = model.get_norm_vals(k_id, norm_params)
        nm_prob, ig_prob = model.event_match_probs((mean, stdv), k_id, norm_params)
        log_prob = np.log(nm_prob*ig_prob)

    print ("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.3f" % (evt_kmer, evt, ref_pos, evt_mean, evt_stdv, lv_mean, lv_stdv, sd_mean, nm_prob, ig_prob, log_prob))
