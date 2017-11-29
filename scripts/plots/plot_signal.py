import h5py
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

if __name__ == "__main__":

    #RAW_ST = 0 + int(sys.argv[2])

    hdf = h5py.File(os.path.abspath(sys.argv[1]),'r')

    raw_reads = hdf['Raw']['Reads']
    event_reads = hdf['Analyses']['EventDetection_000']['Reads']

    
    for r in raw_reads:
        raw = np.array(raw_reads[r]['Signal'])

    for r in event_reads:
        events = event_reads[r]['Events']


    evt_verts = list()
    t = 0
    means = list()
    lengths = list()
    for start, length, mean, stdv in events:
        evt_verts.append((t, mean))
        t += length
        evt_verts.append((t, mean))

        means.append(mean)
        lengths.append(length)

    means = np.array(means)

    scale = means.std() / raw.std()
    shift = means.mean() - (scale * raw.mean())
    raw_norm = raw * scale + shift

    print(np.mean(lengths))

    codes = [Path.MOVETO] + [Path.LINETO]*(len(evt_verts)-1)
    path = Path(evt_verts, codes)
    patch = patches.PathPatch(path, edgecolor='red', facecolor='none', lw=2)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlim(0, t)
    ax.set_ylim(means.min() * 0.95, means.max() * 1.05)

    ax.plot(raw_norm, linewidth=2)
    ax.add_patch(patch)


    plt.show()
