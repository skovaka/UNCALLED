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
from ..dtw.track import Track, RefCoord, ref_coords, method_compare_aln
from ..dtw.read_aln import BcFast5Aln

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
    Opt("track_a", "dotplot", type=str),
    Opt("track_b", "dotplot", nargs="?", type=str),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "track", type=RefCoord),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", type=comma_split),
    Opt(("-C", "--max-chunks"), "read_buffer"),
]

def main(conf):
    """Plot dotplots of alignments from tracks produced by `dtw` or `convert`"""
    matplotlib.use("TkAgg")
    d = Dotplot(conf=conf)
    _RefIndex.index_test(d.tracks[0].index)
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

        if isinstance(self.prms.track_a, Track):
            self.tracks = [self.prms.track_a]
        else:
            self.tracks = [Track(self.prms.track_a, conf=self.conf)]
        #self.colors = [self.prms.styles["color_a"]]

        self.conf.load_config(self.tracks[0].conf)

        if self.conf.dotplot.track_b is not None:
            if isinstance(self.prms.track_b, Track):
                self.tracks.append(self.prms.track_b)
            else:
                self.tracks.append(Track(self.prms.track_b, conf=self.conf))
            #self.colors.append(self.prms.styles["color_b"])

            self.read_ids = self.tracks[0].read_ids & self.tracks[1].read_ids

            if len(self.read_ids) == 0:
                raise RuntimeError("Dotplot tracks must have reads in common")
        else:
            self.read_ids = self.tracks[0].read_ids
        
        self.fast5s = Fast5Reader(reads=self.read_ids, conf=self.tracks[0].conf)

        self.index = self.tracks[0].index
        #self.mm2s = self.tracks[0].mm2s

    def show_all(self):
        for read_id in self.read_ids:
            print(read_id)
            self.show(read_id)

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

    def show(self, read_id=None, cursor=None, fast5_read=None):
        if not self._plot(read_id, cursor, fast5_read):
            return False

        plt.show()
        plt.close()

        return True

    def add_aln(self, aln, focus=False, color="purple", model=None):
        if self.read.id is not None and self.read.id != aln.read_id:
            raise RuntimeError("All Dotplot alignments must be from same read")

        if self.ref_bounds is None:
            self.ref_bounds = aln.ref_bounds

        elif self.ref_bounds.name != aln.ref_bounds.name:
            raise RuntimeError("All Dotplot alignments must be on same reference sequence")
        else:
            self.ref_bounds = RefCoord(
                self.ref_bounds.name,
                min(self.ref_bounds.start, aln.ref_bounds.start), 
                max(self.ref_bounds.end, aln.ref_bounds.end), 
                self.ref_bounds.fwd)

        if focus:
            if len(self.focus) > 0:
                raise RuntimeError("Dotplot can currently only have one focus alignment")
            self.focus.add(len(self.alns))

        self.alns.append((aln, color, model))

    def _plot_signal(self, aln, ax, track):

        samp_min, samp_max = aln.get_samp_bounds()
        sys.stdout.flush()
        raw_norm = self.read.get_norm_signal(samp_min, samp_max)
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

        samps = np.arange(samp_min, samp_max)
        for base, color in enumerate(self.prms.style["base_colors"]):
            ax.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color)

        ax.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, c="black")
        ax.step(aln.aln['start'], model_current, color='white', linewidth=2, where="post")

        skips = aln.aln[np.diff(aln.aln.start, append=-1)==0].index
        ax.vlines(aln.aln.loc[skips, "start"]-1, ymin, ymax, colors="red", linestyles="dashed", linewidth=3)


    def _plot_aln(self, i):
        aln = self.alns[i]
        #if not "refmir" in aln.aln:
        #    aln.calc_refmir()
        #aln.sort_refmir()

        if getattr(aln, "bands", None) is not None:
            self.ax_dot.fill_between(aln.bands['samp'], aln.bands['ref_st']-1, aln.bands['ref_en'], zorder=1, color='#ccffdd', linewidth=1, edgecolor='black', alpha=0.5)

        self.ax_dot.step(aln.aln['start'], aln.aln.index, where="post", linewidth=3,
            **self.prms.style["aln_kws"][i]
        )

        self._plot_signal(aln, self.sig_axs[i], self.tracks[i])


    def set_cursor(self, ref_coord):
        aln,_,_ = self.alns[list(self.focus)[0]]
        refmir = aln.ref_to_refmir(ref_coord)

        i = aln.aln['refmir'].searchsorted(refmir)
        samp = aln.aln.iloc[i]['start'] + aln.aln.iloc[i]['length']/2

        self.cursor = (samp, refmir)

    def _tick_formatter(self, x, pos):
        return self.index.refmir_to_ref(int(x))

    def _load_read(self, read_id, fast5_read=None):
        if read_id is None and fast5_read is None:
            raise ValueError("Must specify read to show")
            #read_id = self.tracks[0].read_aln.read_id

        self.cursor = None

        if fast5_read is None:
            fast5_read = self.fast5s[read_id]
        else:
            read_id = fast5_read.id

        if isinstance(fast5_read, ProcRead):
            self.read = fast5_read
        else:
            self.read = ProcRead(fast5_read, conf=self.conf)

        self.alns = [
            track.load_read(read_id)#, ref_bounds=self.conf.track.ref_bounds)
            for track in self.tracks
        ]

        if np.any([a.empty for a in self.alns]):
            return False

        #self.aln_bc = BcFast5Aln(self.index, self.read, self.mm2s[read_id], refmirs=self.alns[0].refmirs)

        for t in self.tracks:
            t.load_aln_kmers(store=True)
        #self.aln_kmers = [
        #    track.get_aln_kmers(read_id)#, ref_bounds=self.conf.track.ref_bounds)
        #    for track in self.tracks
        #]

        #TODO improve this
        self.ref_bounds = self.alns[0].ref_bounds

        return True
    
    #def _load_layer(self, layer):

    
    @property
    def track_count(self):
        return len(self.tracks)

    def _plot(self, read_id, cursor=None, fast5_read=None):
        if not self._load_read(read_id, fast5_read):
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
            #rax = ax.twinx()
            #rax.set_ylabel(self.tracks[i].name)



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
            samp, refmir = cursor
            self.ax_dot.axvline(samp, **cursor_kw),
            self.ax_dot.axhline(refmir,  **cursor_kw)

        xmin = np.inf
        xmax = -np.inf
        for i in range(self.track_count):
            self._plot_aln(i)
            samp_min, samp_max = self.alns[i].get_samp_bounds()
            xmin = min(xmin, samp_min)
            xmax = max(xmax, samp_max)

        xslop = (xmax-xmin) * 0.01
        self.ax_dot.set_xlim(xmin-xslop, xmax+xslop)

        for ax,layer in zip(self.layer_axs, self.prms.layers):
            for i,aln in enumerate(self.alns):
                ax.step(
                    aln.aln[layer], aln.aln.index+0.5, 
                    where="post",
                    **self.prms.style["aln_kws"][i]
                )

        #if not self.aln_bc.empty:
        #    self.ax_dot.scatter(self.aln_bc.aln['sample'], self.aln_bc.aln.index, color='orange', zorder=2, s=20)

        self.fig.tight_layout()

        return True

