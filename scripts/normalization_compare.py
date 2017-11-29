import sys
import h5py
from scipy.cluster.vq import *
import numpy as np
import os
import scipy.stats as stats


path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+"/..")
from kmer_model import *

model = KmerModel(sys.argv[1])

hdf = h5py.File(sys.argv[2],'r')

bc_events = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']
n = len(bc_events)

means, stdvs, lengths, kmers, ids = [0]*n, [0]*n, [0]*n, [""]*n, [0]*n
#means = [e['mean'] for e in bc_events]
#stdvs = [e['stdv'] for e in bc_events]
#lengths = [e['length'] for e in bc_events]
#kmers = [e['model_state'].astype(str) for e in bc_events]
#ids = [model.kmer_to_i(k) for k in kmers]
for i in range(0, n):
    means[i]   = bc_events[i]['mean']
    stdvs[i]   = bc_events[i]['stdv']
    lengths[i] = bc_events[i]['length']
    kmers[i]   = bc_events[i]['model_state'].astype(str)
    ids[i]     = model.kmer_to_i(kmers[i])

norm_mom = model.get_norm_params(list(zip(means, stdvs, lengths)))

hasT = ['T' in model.i_to_kmer(i) for i in range(0, model.kmer_count)]
has2ndT = [model.i_to_kmer(i)[2] == 'T' for i in range(0, model.kmer_count)]

model_means = [model.get_norm_vals(k)[0] for k in ids]

norm_reg = stats.linregress(model_means, means)[:2]

norm_reg_noT = stats.linregress([m for m, i in zip(model_means, ids) if not hasT[i]], 
                                [m for m, i in zip(means, ids) if not hasT[i]])[:2]

norm_reg_no2ndT = stats.linregress([m for m, i in zip(model_means, ids) if not has2ndT[i]], 
                                   [m for m, i in zip(means, ids) if not has2ndT[i]])[:2]

probs_unnorm = [model.event_match_probs((m, s), k, (1, 0)) for m, s, k in zip(means, stdvs, ids)]

probs_mom = [model.event_match_probs((m, s), k, norm_mom) for m, s, k in zip(means, stdvs, ids)]
probs_reg = [model.event_match_probs((m, s), k, norm_reg) for m, s, k in zip(means, stdvs, ids)]
probs_reg_noT = [model.event_match_probs((m, s), k, norm_reg_noT) for m, s, k in zip(means, stdvs, ids)]
probs_reg_no2ndT = [model.event_match_probs((m, s), k, norm_reg_no2ndT) for m, s, k in zip(means, stdvs, ids)]

print(np.mean([p1*p2 for p1, p2 in probs_unnorm]))
print(np.mean([p1*p2 for p1, p2 in probs_mom]))
print(np.mean([p1*p2 for p1, p2 in probs_reg]))
print(np.mean([p1*p2 for p1, p2 in probs_reg_noT]))
print(np.mean([p1*p2 for p1, p2 in probs_reg_no2ndT]))

print()
print(np.mean([p[0]*p[1] for i, p in zip(ids, probs_mom) if not hasT[i]]))
print(np.mean([p[0]*p[1] for i, p in zip(ids, probs_reg) if not hasT[i]]))
print(np.mean([p[0]*p[1] for i, p in zip(ids, probs_reg_noT) if not hasT[i]] ))
print(np.mean([p[0]*p[1] for i, p in zip(ids, probs_reg_no2ndT) if not hasT[i]]))

