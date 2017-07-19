import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import sys

mean1, stdv1, n = map(float, sys.stdin.readline().strip().split())
var1 = stdv1*stdv1

mean2, stdv2 = map(float, sys.stdin.readline().strip().split())
var2 = stdv2*stdv2

mean_diff = mean2 - mean1
var_diff = var1 + var2
stdv_diff = np.sqrt(var_diff)

for mean, var, stdv in [(mean1, var1, stdv1), (mean2, var2, stdv2), (mean_diff, var_diff, stdv_diff)]:
    x = np.linspace(mean-3*var,mean+3*var, 100)
    plt.plot(x, mlab.normpdf(x, mean, stdv))
    plt.scatter([mean - stdv, mean + stdv], [0, 0])
plt.scatter([mean_diff - 2*stdv_diff, mean_diff + 2*stdv_diff], [0, 0])

print (stats.ttest_ind_from_stats(mean_diff, stdv_diff, n, 0, 0, 2, equal_var = False))
print (stats.ttest_ind_from_stats(mean1, stdv1, n, mean2, stdv2, 1000, equal_var = False))
print (mean_diff, stdv_diff)
print (mean_diff / np.sqrt(var_diff / n))

plt.show()


