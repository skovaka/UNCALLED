import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
import re
from glob import glob

param_names = {'m' : "model",
               'r' : "reference",
               'l' : "read_list",
               'o' : "out_prefix",
               't' : "tally_dist",
               's' : "min_seed_len",
               'a' : "anchor_len",
               'i' : "max_ignores",
               'k' : "max_skips",
               'y' : "stay_frac",
               'E' : "extend_evpr",
               'A' : "anchor_evpr",
               'S' : "seed_pr",
               'Y' : "stay_pr"}


param_re = re.compile(r"(-?\d+(\.\d+)?)(\.[A-Za-z]|_)")

pattern = sys.argv[1]

y_type = sys.argv[2]

if y_type == "time":
    y_time = True
    y_frac = y_ratio = y_first = False

elif y_type == "frac":
    y_frac = True
    y_time = y_ratio = y_first = False

elif y_type == "ratio":
    y_ratio = True
    y_time = y_frac = y_first = False

elif y_type == "firstseed":
    y_first = True
    y_time = y_frac = y_ratio = False

else:
    print ("Unrecognized y type %s" % y_type)
    sys.exit()

filenames = glob(pattern)
print(filenames)
i = pattern.find("*")
p = pattern[i-1]
x_axis = param_names[p]

data = dict()

for fname in filenames:
    print(fname[i:])
    m = param_re.search(fname[i:])
    x = float(m.group(1))

    y = list()

    infile = open(fname)

    while True:
        infile.readline()
        line = infile.readline()

        tabs = line.split()
        if len(tabs) == 0:
            break
        event_ct = int(tabs[-3])

        top_fwd = infile.readline().strip().split()
        top_rev = None
        next_fwd = None
        next_rev = None

        eof = True

        if len(top_fwd) > 0 and top_fwd[0][0] != "=":

            for line in infile:
                tabs = line.strip().split()
                eof = False
                if tabs[0] == "fwd":
                    if next_fwd == None and float(tabs[-1]) < 0.9:
                        next_fwd = tabs
                else:
                    if tabs[0] == "rev":
                        top_rev = tabs
                    break

            if eof:
                break
            
            if top_rev != None:
                for line in infile:
                    tabs = line.strip().split()
                    if tabs[0] == "rev":
                        if next_rev == None and float(tabs[-1]) < 0.9:
                            next_rev = tabs
                    else:
                        break

        print (line)

        if top_fwd != None and top_rev != None:
            if int(top_fwd[1]) > int(top_rev[1]):
                best = top_fwd
                next_best = next_fwd
            else:
                best = top_rev
                next_best = next_rev
        elif top_fwd:
            best = top_fwd
            next_best = next_fwd
        elif top_rev:
            best = top_rev
            next_best = next_rev
        else:
            best = None

        if best != None:
            time = float(tabs[1])
            
            if y_frac:
                length = int(best[1])
                st, en = map(int, best[2].split("-"))
                y.append(length / (en - st + 1))

            elif y_time:
                y.append(time / event_ct)

            elif y_ratio:
                if next_best == None:
                    y.append(-1)
                else:
                    y.append(float(next_best[-2]))

            elif y_first:
                st, en = map(int, best[3].split("-"))
                y.append(event_ct - en)

    for j in range(0, len(y)):
        if y[j] < 0:
            y[j] = max(y)
    
    y_mean = np.mean(y)
    y_stdv = np.std(y)
    #data[x] = (y_mean, y_stdv)
    data[x] = y

    print("%f\t%f" % (x, y_mean))

xs = list()
y_means = list()
y_stdvs = list()

d = list()
for x in sorted(data.keys()):
    d.append(data[x])
#    y_mean, y_stdv = data[x]
#    xs.append(x)
#    y_means.append(y_mean)
#    y_stdvs.append(y_stdv)

#slope, incp, r, p, ster = stats.linregress(xs, ys)

#print ("slope=%.2f, intercept=%.2f, p=%.4f, r=%.4f" % (slope, incp, p, r))

#plt.errorbar(xs, y_means, yerr=y_stdvs, capsize=5)
plt.boxplot(d)

plt.xlabel(x_axis)
if y_frac:
    plt.ylabel("Fraction read aligned")

elif y_time:
    plt.ylabel("Time / Event (sec)")

elif y_ratio:
    plt.ylabel("Best / Next Best ratio")


plt.show()
