import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats

probs = list()
tp_fracs = list()
fp_fracs = list()
tp_fp_rat = list()

tp_in = open(sys.argv[1])
fp_in = open(sys.argv[2])

for line in tp_in:
    prob, frac = map(float, line.strip().split())
    probs.append(prob)
    tp_fracs.append(frac)

for line in fp_in:
    prob, frac = map(float, line.strip().split())
    fp_fracs.append(frac)

for tp, fp in zip(tp_fracs, fp_fracs):
    tp_fp_rat.append(fp - tp)

plt.plot(probs, fp_fracs)
plt.plot(probs, tp_fracs)
plt.plot(probs, tp_fp_rat)

#plt.title("TP and FP seeds removed at ")
plt.xlabel("Log min event probability")
plt.ylabel("Proportion of seeds removed")
plt.legend(['fp', 'tp', 'diff'])
plt.show()

