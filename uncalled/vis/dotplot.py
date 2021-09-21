import time
import os 
import sys
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter
from matplotlib.patches import Rectangle
import types

from .. import nt, config

from ..sigproc import ProcRead
from ..fast5 import Fast5Reader, parse_read_ids
from ..pafstats import parse_paf
from ..dtw.track import AlnTrack, RefCoord, ref_coords, method_compare_aln
from ..dtw.track_io import TrackIO

from _uncalled import _RefIndex

class DotplotParams(config.ParamGroup):
    _name = "dotplot"
DotplotParams._def_params(
    ("track_a", None, None, "DTW aligment track containing reads to plot"),
    ("track_b", None, None, "DTW aligment track containing reads to plot"),
    ("layers", ["current", "length"], list, "Layers to plot in side panels"),
    ("style", {
        "rc" : {
            "figure.figsize" : (15, 10), 
            "axes.labelsize" : 14, 
            "axes.titlesize" : 16},
        "aln_kws" : [
            {"color" : "purple", "alpha" : 1},
            {"color" : "royalblue", "alpha" : 0.75}],
        "bc_color": "orange",
        "base_colors" : ["#80ff80", "#8080ff", "#ffbd00", "#ff8080"], #A,C,G,T
    }, dict, "Plotting style options")
)
#Config._EXTRA_GROUPS["dotplot"] = DotplotParams #TODO put in ParamGroup con

def comma_split(s):
    return s.split(",")

from ..argparse import Opt
OPTS = [
    #Opt("track_a", "dotplot", type=str),
    #Opt("track_b", "dotplot", nargs="?", type=str),
    Opt("input", "track_io", nargs="+"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "track_io", type=RefCoord),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", type=comma_split),
    Opt(("-C", "--max-chunks"), "read_buffer"),
]

def main(conf):
    """Plot dotplots of individual DTW alignments"""
    matplotlib.use("TkAgg")
    d = Dotplot(conf=conf)
    if conf.out_prefix is None:
        d.show_all()
    else:
        d.save_all(conf.out_prefix, conf.out_format)


class Dotplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("dotplot", *args, **kwargs)

        #matplotlib.use("TkAgg")
        matplotlib.rcdefaults()
        plt.style.use(['seaborn'])
        matplotlib.rcParams.update(self.prms.style["rc"])

        self.conf.track.load_mat = False
        self.conf.track.layers = self.conf.dotplot.layers

        self.track_io = TrackIO(conf=self.conf)
        self.fast5s = self.track_io.get_fast5_reader()#Fast5Reader(reads=self.read_ids, index=fast5_index, conf=self.tracks[0].conf)
        self.conf.load_config(self.track_io.conf)

        self.index = self.track_io.index
        #self.mm2s = self.tracks[0].mm2s

    def show_all(self):
        for read_id, tracks in self.track_io.iter_reads():
            self.read_id = read_id
            self.tracks = tracks
            print("TRACKSSS", tracks)
            self.show(tracks)

    def save_all(self, out_prefix, fmt):
        for read_id in self.read_ids:
            if not self._plot(read_id):
                continue
            print(read_id)

            suffix = read_id + "." + fmt
            if os.path.isdir(out_prefix):
                out_fname = os.path.join(out_prefix, suffix)
            else:
                out_fname = out_prefix + suffix

            self.fig.set_size_inches(8, 6)
            self.fig.savefig(out_fname, dpi=200)
            plt.close()

    def show(self, tracks=None, cursor=None, fast5_read=None):
        if not self._plot(cursor, fast5_read):
            return False

        plt.show()
        plt.close()

        return True

    def _plot_signal(self, ax, track):

        dtw = track.layers["dtw"]

        samp_min = dtw["start"].min()
        max_i = dtw["start"].argmax()
        samp_max = dtw["start"].iloc[max_i] + dtw["length"].iloc[max_i]
        #samp_min, samp_max = aln.get_samp_bounds()

        raw_norm = self.read.get_norm_signal(samp_min, samp_max)

        kmers = track.coords.kmers[dtw.index.get_level_values(0)]
        model_current = self.track_io.model[kmers]

        ymin = min(np.min(model_current), np.min(raw_norm[raw_norm>0]))
        ymax = max(np.max(model_current), np.max(raw_norm[raw_norm>0]))

        starts = dtw['start']
        ends = starts+dtw['length']

        aln_bases = nt.kmer_base(kmers, 2)
        samp_bases = np.full(samp_max-samp_min, -1)

        for i in range(len(dtw)):
            st = int(dtw.iloc[i]['start'] - samp_min)
            en = int(st + dtw.iloc[i]['length']) - 1
            samp_bases[st:en] = aln_bases[i]

        samps = np.arange(samp_min, samp_max)
        for base, color in enumerate(self.prms.style["base_colors"]):
            ax.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color)

        ax.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, c="black")
        ax.step(dtw['start'], model_current, color='white', linewidth=2, where="post")

        return samp_min, samp_max

        #skips = aln.dtw[np.diff(aln.dtw.start, append=-1)==0].index
        #ax.vlines(aln.dtw.loc[skips, "start"]-1, ymin, ymax, colors="red", linestyles="dashed", linewidth=3)


    def _plot_aln(self, i):
        track = self.tracks[i]
        #aln = self.alns[i]
        #if not "mref" in aln.dtw:
        #    aln.calc_mref()
        #aln.sort_mref()

        #if getattr(aln, "bands", None) is not None:
        #    self.ax_dot.fill_between(aln.bands['samp'], aln.bands['ref_st']-1, aln.bands['ref_en'], zorder=1, color='#ccffdd', linewidth=1, edgecolor='black', alpha=0.5)

        mrefs = track.layers.index.get_level_values("mref")

        self.ax_dot.step(track.layers["dtw"]["start"], mrefs, where="post", linewidth=3,
            **self.prms.style["aln_kws"][i]
        )

        return self._plot_signal(self.sig_axs[i], self.tracks[i])


    def set_cursor(self, ref_coord):
        aln,_,_ = self.alns[list(self.focus)[0]]
        mref = aln.ref_to_mref(ref_coord)

        i = aln.dtw['mref'].searchsorted(mref)
        samp = aln.dtw.iloc[i]['start'] + aln.dtw.iloc[i]['length']/2

        self.cursor = (samp, mref)

    def _tick_formatter(self, x, pos):
        return self.index.mref_to_ref(int(x))

    def _load_read(self):#, read_id, fast5_read=None):
        #if read_id is None and fast5_read is None:
        #    raise ValueError("Must specify read to show")
        #    #read_id = self.tracks[0].read_aln.read_id

        self.cursor = None

        fast5_read = self.fast5s[self.read_id]

        if isinstance(fast5_read, ProcRead):
            self.read = fast5_read
        else:
            print("POOR MODEL", self.conf.pore_model.name)
            self.read = ProcRead(fast5_read, conf=self.conf)

        empty = False
        for track in self.tracks:
            empty = empty or len(track.layers) == 0

        if empty:
            return False

        print("TRACK", self.tracks)

        #TODO improve this
        self.ref_bounds = self.tracks[0].aln_ref_coord(self.tracks[0].alignments.index[0])

        return True
    
    #def _load_layer(self, layer):

    
    @property
    def track_count(self):
        return len(self.tracks)

    def _plot(self, cursor=None, fast5_read=None):
        if not self._load_read():#read_id, fast5_read):
            return False

        widths=[5]
        for _ in self.prms.layers:
            widths.append(1)

        heights=[1] * self.track_count
        heights.append(3)

        self.fig = plt.figure()
        gspec = self.fig.add_gridspec(
                ncols=len(widths), nrows=len(heights), 
                width_ratios=widths,
                height_ratios=heights 
        )

        track_count = len(self.tracks)

        self.sig_axs = [self.fig.add_subplot(gspec[i,0]) 
                        for i in range(track_count)]

        self.ax_dot = self.fig.add_subplot(gspec[track_count,0])

        self.layer_axs = [self.fig.add_subplot(gspec[track_count,l+1])
                          for l in range(len(self.prms.layers))]

        self.sig_axs[0].xaxis.tick_top()
        self.sig_axs[0].xaxis.set_label_position("top")
        for ax in self.sig_axs[1:]:
            ax.xaxis.set_major_formatter(NullFormatter())

        for i,ax in enumerate(self.sig_axs):
            ax.get_shared_x_axes().join(self.ax_dot, ax)
            ax.set_ylabel(
                "Current (pA)", 
                color=self.prms.style["aln_kws"][i]["color"])

        self.sig_axs[0].set_title(self.read.id)
        self.ax_dot.set_xlabel("Raw Sample")

        for i,ax in enumerate(self.layer_axs):
            ax.get_shared_y_axes().join(self.ax_dot, ax)
            ax.yaxis.set_major_formatter(NullFormatter())
            ax.set_xlabel(self.prms.layers[i])

        self.ax_dot.set_ylabel("%s (%s)" % (self.ref_bounds.name, "+" if self.ref_bounds.fwd else "-"))

        self.ax_dot.yaxis.set_major_formatter(FuncFormatter(self._tick_formatter))
        
        if cursor is not None:
            cursor_kw = {
                'color' : 'red', 
                'alpha' : 0.5
            }
            samp, mref = cursor
            self.ax_dot.axvline(samp, **cursor_kw),
            self.ax_dot.axhline(mref,  **cursor_kw)

        xmin = np.inf
        xmax = -np.inf
        for i in range(self.track_count):
            samp_min, samp_max = self._plot_aln(i)
            xmin = min(xmin, samp_min)
            xmax = max(xmax, samp_max)

        xslop = (xmax-xmin) * 0.01
        self.ax_dot.set_xlim(xmin-xslop, xmax+xslop)

        for ax,layer in zip(self.layer_axs, self.prms.layers):
            for i,track in enumerate(self.tracks):
                mrefs = track.layers.index.get_level_values("mref")
                ax.step(
                    track.layers["dtw"][layer], mrefs+0.5, 
                    where="post",
                    **self.prms.style["aln_kws"][i]
                )

        #if not self.aln_bc.empty:
        #    self.ax_dot.scatter(self.aln_bc.aln['sample'], self.aln_bc.aln.index, color='orange', zorder=2, s=20)

        self.fig.tight_layout()

        return True

