import sys
import h5py
import matplotlib.pyplot as plt
from scipy.cluster.vq import *
import numpy as np
import os
import scipy.stats as stats


path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+"/..")
from kmer_model import *

model = KmerModel(sys.argv[1])


hdf = h5py.File(sys.argv[2],'r')

events = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']

l = int(len(events)/3)
events = events[l:2*l]

events = [(e['mean'], e['stdv'], e['model_state'].astype(str)) for e in events]

#events = list(filter(lambda e: "T" not in e[2], events))


norm = model.get_norm_params(events)

event_means = [e[0] for e in events]
event_kmers = [model.kmer_to_i(e[2]) for e in events]
model_means = [model.get_norm_vals(k)[0] for k in event_kmers]

slope, incp, r, p, ster = stats.linregress(model_means, event_means)
reg_norm = (slope, incp)
print(norm, reg_norm)

probs1 = [model.event_match_probs((m, s), model.kmer_to_i(k), norm) for m, s, k in events]
probs2 = [model.event_match_probs((m, s), model.kmer_to_i(k), (1, 0)) for m, s, k in events]
probs3 = [model.event_match_probs((m, s), model.kmer_to_i(k), reg_norm) for m, s, k in events]

print(np.mean([p1*p2 for p1, p2 in probs1]))
print(np.mean([p1*p2 for p1, p2 in probs2]))
print(np.mean([p1*p2 for p1, p2 in probs3]))

colors = ['blue' if e[2][1] == 'T' else 'blue' for e in events]

plt.axis([0, 160, 0, 160])


plt.scatter(event_means, model_means, alpha=0.1, c=colors)
#plt.plot([0, 160], [incp, incp+160*slope], c='red')
plt.xlabel("Read Event Mean (pA)")
plt.ylabel("Model Event Mean (pA)")
plt.show()


