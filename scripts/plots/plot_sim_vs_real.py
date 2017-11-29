
import sys, os
path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+"/..")

import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from match_sim import *

model = KmerModel(sys.argv[1])
seq_len = int(sys.argv[2])
N = int(sys.argv[3])

event_thresh = -5.29
seed_thresh  = -3.35

event_probs = np.arange(-9.1, -3.5, 0.1)

#for seq_len in range(13, 70, 8):
seq = "".join([np.random.choice(["A", "C", "G", "T"]) for _ in range(seq_len)])
print(seq)
data = list()

d = 0

data.append([0] * len(event_probs))
for i in range(0, N):
    seq2 = mutate_seq(seq, d)
    seed_prob, min_prob = match_prob(model, seq, seq2)
    
    for p in range(0, len(event_probs)):
         data[-1][p] += min_prob >= event_probs[p]

data[-1] = [float(d) / N for d in data[-1]]

for d in data:
    plt.plot(event_probs, d)

plt.xlabel("Log min event probability")
plt.ylabel("Proportion of matching sequences")
#plt.legend(["Edit distance %d" % d for d in range(0, len(data))])

tp_fracs = list()
tp_probs = list()

tp_in = open(sys.argv[4])

for line in tp_in:
    prob, frac = map(float, line.strip().split())
    tp_probs.append(prob)
    tp_fracs.append(1-frac)

plt.plot(tp_probs, tp_fracs)

#plt.title("TP and FP seeds removed at ")
#plt.xlabel("Log min event probability")
#plt.ylabel("Proportion of seeds removed")
#plt.legend(['fp', 'tp', 'diff'])
plt.show()

