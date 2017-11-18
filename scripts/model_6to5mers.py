from kmer_model import *
import sys
import numpy as np
import matplotlib.pyplot as plt

model_6mers = KmerModel(sys.argv[1])
model_5mers = KmerModel(k=5)

lv_diffs = list() 
sd_diffs = list()

for i in range(0, model_5mers.kmer_count):
    kmer5 = model_5mers.i_to_kmer(i)

    ids6 = [model_6mers.kmer_to_i(b+kmer5) for b in model_5mers.BASES]

    vals = [model_6mers.get_norm_vals(k) for k in ids6]
    lvmean5 = np.mean([lvmean for lvmean, lvstdv, sdmean, sdstdv in vals])
    lvstdv5 = np.mean([lvstdv for lvmean, lvstdv, sdmean, sdstdv in vals])
    sdmean5 = np.mean([sdmean for lvmean, lvstdv, sdmean, sdstdv in vals])
    sdstdv5 = np.mean([sdstdv for lvmean, lvstdv, sdmean, sdstdv in vals])

    print ("%s\t%.6f\t%.6f\t%.6f\t%.6f\t0" % (kmer5, lvmean5, lvstdv5, sdmean5, sdstdv5))

    for lvmean6, lvstdv6, sdmean6, sdstdv6 in vals:
        lv_diffs.append(abs(lvmean5 - lvmean6) / lvstdv6) 
        sd_diffs.append(abs(sdmean5 - sdmean6) / sdstdv6)

f = len(list(filter(lambda d: d < 1, lv_diffs))) / float(len(lv_diffs))

sys.stderr.write("Fraction of diffs < 1 stdv: %f\n" % f)

plt.hist(lv_diffs, bins=50)
plt.xlabel("Number of stdvs each 6-mer mean is from it's 5-mer mean", fontsize=25)
plt.show()
