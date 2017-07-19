import sys
import matplotlib.pyplot as plt
from scipy.cluster.vq import *
import numpy as np

sys.stdin.readline()

k = int(sys.argv[1])

COLORS = ['blue', 'red', 'green', 'black', 'orange', 'magenta', 'cyan', 'yellow'] 

kmers = list()
lv_means = list()
lv_stdvs = list()
sd_means = list()
sd_stdvs = list()
ig_lambdas = list()
weights = list()

obs = list()

for line in sys.stdin:
    tabs = line.strip().split()
    kmer = tabs[0]
    lv_mean, lv_stdv, sd_mean, sd_stdv, ig_lambda, weight = map(float, tabs[1:])
    
    kmers.append(kmer)
    lv_means.append(lv_mean)
    lv_stdvs.append(lv_stdv)
    sd_means.append(sd_mean)
    sd_stdvs.append(sd_stdv)
    ig_lambdas.append(ig_lambda)
    weights.append(weight)

    #obs.append([lv_mean, lv_stdv, sd_mean, sd_stdv])
    obs.append([lv_mean, sd_mean])
    #obs.append([sd_mean, sd_stdv])

obs = whiten(np.array(obs))

centroids, distortion = kmeans(obs, k)
codes, dists = vq(obs, centroids)

base_counts = list()
for i in range(0, k):
    base_counts.append([0, 0, 0, 0])

for i in range(0, len(kmers)):
    base_counts[codes[i]][0] += kmers[i].count('A')
    base_counts[codes[i]][1] += kmers[i].count('C')
    base_counts[codes[i]][2] += kmers[i].count('G')
    base_counts[codes[i]][3] += kmers[i].count('T')


print ("\tA\tC\tG\tT")
for i in range(0, k):
    percents = tuple(map(lambda c: 100.0*c/sum(base_counts[i]),  base_counts[i]))
    print ("%s\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%" % ((COLORS[i],)+percents))

#plt.scatter(sd_stdvs, sd_means, c=list(map(lambda c: COLORS[c], codes)))
plt.scatter(sd_means, lv_means, c=list(map(lambda c: COLORS[c], codes)))
plt.show()










