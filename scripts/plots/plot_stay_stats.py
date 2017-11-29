import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

#infile = open(sys.argv[1])
infile = sys.stdin
target = int(sys.argv[1])

seed_lens = map(int, infile.readline().strip().split())

data = [list() for _ in seed_lens]

readnm = -1

for line in infile:
    if line[0] == "=":
        readnm += 1
        continue
    
    if readnm < target:
        continue
    
    if readnm > target:
        break
    
    s = 0
    for f in line.strip().split():
        if f != "NA":
            data[s].append(float(f))
        s += 1

xs = np.arange(0, 1, 0.001)

for fracs in data:
    fracs = sorted(fracs)
    
    ys = list()

    count = 0
    i = 0
    above = False
    for x in xs:
        while i < len(fracs) and fracs[i] < x:
            count += 1
            i += 1

        ys.append(count / len(fracs))

        if not above and ys[-1] > 0.99:
            above = True
            print(x)

    plt.plot(xs, ys)

plt.legend(['8', '64'])
plt.show()



#plt.plot(tp_probs, tp_fracs)

#plt.title("TP and FP seeds removed at ")
#plt.xlabel("Log min event probability")
#plt.ylabel("Proportion of seeds removed")
#plt.legend(['fp', 'tp', 'diff'])
#plt.show()

