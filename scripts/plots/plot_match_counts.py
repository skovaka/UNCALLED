import sys
import matplotlib.pyplot as plt
from scipy.cluster.vq import *
import numpy as np
import os
import h5py

path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+"/..")
from kmer_model import *

model = KmerModel(sys.argv[1])
d = sys.argv[2]

i = 0

hdf = h5py.File(os.path.abspath(d),'r')

reads = hdf['Analyses']['EventDetection_000']['Reads']

counts = list()

for r in reads:
    events = [ (e[2], e[3], e[1]) for e in reads[r]['Events'] ]
    n = model.get_norm_params(events)

    for e in events:
        ct = 0
        for k in range(0, model.kmer_count):
            ng,ig = model.event_match_probs(e[:2], k, n)
            l = np.log(ig*ng)
            if l > -3.75:
                ct += 1
        print (ct)
        sys.stdout.flush()

#plt.histogram(counts)

#plt.imshow(heatmap.T, extent=extent, origin='lower')
#plt.scatter(lv_means, sd_means, s=lv_stdvs, c=(0, 0, 1), alpha=0.2, linewidths=0)
#plt.xlabel("Expected Mean (pA)", fontsize=15)
#plt.ylabel("Expected Stdv (pA)", fontsize=15)
#plt.title('Pore Model Event Values', fontsize=15)
#plt.show()










