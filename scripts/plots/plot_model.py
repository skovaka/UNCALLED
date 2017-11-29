import sys
import matplotlib.pyplot as plt
from scipy.cluster.vq import *
import numpy as np
import os

path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+"/..")
from kmer_model import *


model = KmerModel(sys.argv[1])

lv_means = list()
sd_means = list()
lv_stdvs = list()
for i in range(0, model.kmer_count):
    lv_mean, lv_stdv, sd_mean, sd_stdv = model.get_norm_vals(i)

    lv_means.append(lv_mean)
    sd_means.append(sd_mean)
    lv_stdvs.append((lv_stdv+sd_stdv)*20)

heatmap, xedges, yedges = np.histogram2d(lv_means, sd_means, weights=lv_stdvs, bins=100, normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

#plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.scatter(lv_means, sd_means, s=lv_stdvs, c=(0, 0, 1), alpha=0.2, linewidths=0)
plt.xlabel("Expected Mean (pA)", fontsize=25)
plt.ylabel("Expected Stdv (pA)", fontsize=25)
plt.title('Pore Model Event Values', fontsize=25)
plt.show()










