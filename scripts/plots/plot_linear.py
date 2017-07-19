import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats

xs = list()
ys = list()
for line in sys.stdin:
    x, y = line.strip().split()
    xs.append(float(x))
    ys.append(float(y))

slope, incp, r, p, ster = stats.linregress(xs, ys)

print ("slope=%.2f, intercept=%.2f, p=%.4f, r=%.4f" % (slope, incp, p, r))

plt.scatter(xs, ys)

plt.show()

