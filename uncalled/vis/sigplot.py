import time
import os 
import sys
from collections import defaultdict
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter
from matplotlib.patches import Rectangle
import types

from .. import nt, config

from ..dtw.dtw import Fast5Processor
from ..fast5 import Fast5Reader, parse_read_ids
from ..pafstats import parse_paf
from ..dtw import RefCoord, Track, ref_coords, BcFast5Aln
from ..dtw.track import  method_compare_aln

from _uncalled import _RefIndex

class SigplotParams(config.ParamGroup):
    _name = "sigplot"
SigplotParams._def_params(
    ("tracks", list, None, "DTW aligment tracks"),
    ("reads", list, None, "Reads to plot"),
    ("connectors", bool, True, "Reads to plot"),
    #("axes", None, None, "Axes to plot each read"),
    ("style", {
        "rc" : {
            "figure.figsize" : (15, 10), 
            "axes.labelsize" : 14, 
            "axes.titlesize" : 16},
        "base_colors" : ["#80ff80", "#8080ff", "#ffbd00", "#ff8080"], #A,C,G,T
    }, dict, "Plotting style options")
)

def comma_split(s):
    return s.split(",")

from ..argparse import Opt
OPTS = [
    Opt("tracks", "sigplot", type=str, nargs="+"),
    Opt(("-l", "--reads"), "sigplot", type=parse_read_ids),
    Opt(("-R", "--ref-bounds"), "track", type=RefCoord),
    Opt(("-C", "--max-chunks"), "read_buffer"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
]

def main(conf):
    """Plot dotplots of alignments from tracks produced by `dtw` or `convert`"""
    Sigplot(conf=conf).show()
    #if conf.out_prefix is None:
    #    d.show_all()
    #else:
    #    d.save_all(conf.out_prefix, conf.out_format)

class Sigplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("sigplot", *args, **kwargs)

        matplotlib.use("TkAgg")
        matplotlib.rcdefaults()
        plt.style.use(['seaborn'])
        matplotlib.rcParams.update(self.prms.style["rc"])

        self.conf.track.load_mat = False
        self.conf.track.layers = self.conf.dotplot.layers

        self.tracks = list()
        for t in self.prms.tracks:
            if isinstance(t, str):
                t = Track(t, conf=self.conf)
            t.fast5s = Fast5Processor(conf=t.conf)
            self.tracks.append(t)

        self.conf.load_config(self.tracks[0].conf)

        if self.prms.reads is None:
            self.read_ids = self.tracks[0].read_ids
            for t in self.tracks:
                self.read_ids = self.read_ids & t.read_ids
        else:
            self.read_ids = self.prms.reads

        self.index = self.tracks[0].index

    def show(self, read_ids=None):
        if read_ids is not None:
            self.read_ids = read_ids

        if not self._plot():
            return False

        plt.show()
        plt.close()

        return True

    def _load_reads(self):
        self.reads = dict()
        self.alns = list()
        self.aln_count = 0

        for read_id in self.read_ids:
            for track in self.tracks:
                self._load_aln(track, read_id)
        return True

    def _load_aln(self, track, read_id):
        if not read_id in track: return
        aln = track.load_read(read_id)
        if aln.empty: return

        self.alns.append(aln)
        self.aln_count += 1

        if read_id not in self.reads:
            self.reads[read_id] = (track, track.fast5s[read_id])

        return True

    @property
    def track_count(self):
        return len(self.tracks)

    def _plot(self):
        if not self._load_reads():
            return False

        self.fig, self.axs = plt.subplots(len(self.alns), 1)
        #for ax in self.axs[1:]:
        #    ax.get_shared_x_axes().join(self.axs[0], ax)

        for i in range(len(self.alns)):
            self.plot_aln(i)
        return True

    def plot_aln(self, i):
        ax = self.axs[i]
        aln = self.alns[i]
        read_id = aln.read_id

        ax.xaxis.set_major_formatter(NullFormatter())
        #ax.set_title(read_id)
        #ax.set_xlabel("Raw Sample")
        #ax.set_ylabel("Current (pA)")

        samp_min, samp_max = aln.get_samp_bounds()

        xshift = samp_min + (samp_max - samp_min) / 2
        #yshift = 

        track, read = self.reads[read_id]
        raw_norm = read.get_norm_signal(samp_min, samp_max)
        model_current = track.model[aln.aln['kmer']]

        ymin = min(np.min(model_current), np.min(raw_norm[raw_norm>0]))
        ymax = max(np.max(model_current), np.max(raw_norm[raw_norm>0]))

        starts = aln.aln['start']
        ends = starts+aln.aln['length']

        aln_bases = nt.kmer_base(aln.aln['kmer'], 2)
        samp_bases = np.full(samp_max-samp_min, -1)

        for i in range(len(aln.aln)):
            st = int(aln.aln.iloc[i]['start'] - samp_min)
            en = int(st + aln.aln.iloc[i]['length']) - 1
            samp_bases[st:en] = aln_bases[i]

        samps = np.arange(samp_min, samp_max) - xshift
        for base, color in enumerate(self.prms.style["base_colors"]):
            ax.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color)

        ax.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, c="black")

        ax.step(aln.aln['start']-xshift, model_current, color='white', linewidth=2, where="post")

        skips = aln.aln[np.diff(aln.aln.start, append=-1)==0].index
        ax.vlines(aln.aln.loc[skips, "start"]-xshift-1, ymin, ymax, colors="red", linestyles="dashed", linewidth=3)

        self.fig.tight_layout()

        return True

