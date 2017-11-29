import sys, os
path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+"/..")

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from match_sim import *

model = KmerModel(sys.argv[1])
seq_len = int(sys.argv[2])
N = int(sys.argv[3])

event_thresh = -5.29
seed_thresh  = -3.35

#seq = filter(lambda s: s[0] != ">", fasta.readlines())
#seq = "".join(map(lambda s: s.strip(), seq))

if len(sys.argv) > 4:
    dist_max = int(sys.argv[4])
else:
    dist_max = 3

event_probs = np.arange(-11, -3.5, 0.1)

#for seq_len in range(13, 70, 8):
seq = "".join([np.random.choice(["A", "C", "G", "T"]) for _ in range(seq_len)])
print(seq)
data = list()

for d in range(0, dist_max+1):
    data.append([0] * len(event_probs))
    for i in range(0, N):
        seq2 = mutate_seq(seq, d)
        seed_prob, min_prob = match_prob(model, seq, seq2)
        
        for p in range(0, len(event_probs)):
             data[-1][p] += min_prob >= event_probs[p]

    data[-1] = [float(d) / N for d in data[-1]]

plt.clf()

for d in data:
    plt.plot(event_probs, d)


plt.title("Seed len %dbp" % (seq_len))
plt.xlabel("Log min event probability")
plt.ylabel("Proportion of matching sequences")
plt.legend(["Edit distance %d" % d for d in range(0, len(data))])

plt.show()
#plt.axis([event_probs[0], event_probs[-1], 0, 1])
#plt.savefig("simulated_match_dist/%d_%d.png" % (seq_len, N), dpi=200)
