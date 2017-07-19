import numpy as np
import matplotlib.pyplot as plt
import sys


ref_file = open(sys.argv[1])
ref_locs = list()
ref_means = list()
ref_stdvs = list()

for line in ref_file:
    i, kmer, mean, stdv = line.strip().split()
    ref_locs.append(int(i))
    ref_means.append(float(mean))
    ref_stdvs.append(float(stdv))

scale = 0.915338
shift = 11.0352

read_file = open(sys.argv[2])
read_locs = list()
read_means = list()
read_stdvs = list()
read_file.readline()
for line in read_file:
    i, mean, stdv = line.strip().split()
    read_locs.append(float(i)*0.444399839)
    read_means.append(float(mean)-shift)
    read_stdvs.append(float(stdv))

plt.plot(ref_locs, ref_means)
plt.plot(read_locs, read_means)

plt.show()

