import h5py
import sys
import numpy as np
import argparse
import os

if __name__ == "__main__":
    hdf = h5py.File(os.path.abspath(sys.argv[1]),'r')

    path = list()

    while True:

        h = hdf
        for a in path:
            h = h[a]

        print ("-------------")
        if hasattr(h, "keys"):
            i = 0
            for k in h.keys():
                print ("%02d %s" % (i, k))
                i += 1
        else:
            for x in h[:20]:
                print (x)
        print ("-------------")

        sys.stdout.write("> ")
        line = sys.stdin.readline()

        if len(line.strip()) == 0:
            if len(path) == 0:
                break

            path.pop()
            continue

        try:
            i = int(line.strip())

            if i < 0 or i >= len(h.keys()):
                print ("Nope")
                continue
            else:
                print(list(h.keys()), i)
                path.append(list(h.keys())[i])
        except:
            print ("Not right")

    print ()
