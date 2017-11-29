import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats

d1 = list()
for line in open(sys.argv[1]):
    d1.append(int(line.strip()))

d2 = list()
for line in open(sys.argv[2]):
    d2.append(int(line.strip()))

# the histogram of the data
plt.hist(d2, 40, facecolor='red', alpha=0.75, range=(0, 2000))
plt.hist(d1, 40, facecolor='blue', alpha=0.75, range=(0, 2000))
#plt.xticks(np.arange(0, max(data), 500))
#plt.axis([40, 160, 0, 0.1])
plt.xlabel("Number of Matched K-Mers", fontsize=30)
plt.ylabel("Event Count", fontsize=30)
plt.title('K-Mers Matched per Event', fontsize=30)
#perm (min prob. 4.5*10^-5)
#strict (min prob 2.4*10^-2)
plt.legend(["Anchor Threshold", "Extend Threshold"], fontsize=30)

# add a 'best fit' line

plt.show()

