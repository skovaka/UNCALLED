import time
import numpy as np
import os 

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter

from .. import nt, BwaIndex, config

from ..sigproc import ProcRead
from ..fast5 import Fast5Reader, parse_read_ids
from ..pafstats import parse_paf
from . import RefCoord, Track, ref_coords, BcFast5Aln, method_compare_aln

class DotplotParams(config.ParamGroup):
    _name = "dotplot"
DotplotParams._def_params(
    ("track_a", None, str, "DTW aligment track containing reads to plot"),
    ("track_b", None, str, "DTW aligment track containing reads to plot"),
)
#Config._EXTRA_GROUPS["dotplot"] = DotplotParams #TODO put in ParamGroup con

from ..config import Opt
OPTS = [
    Opt("track_a", "dotplot"),
    Opt("track_b", "dotplot", nargs="?"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "track", type=RefCoord),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-C", "--max-chunks"), "read_buffer"),
]

def main(conf):
    """Plot dotplots of alignments from tracks produced by `align` or `convert`"""
    Dotplot(conf=conf).show_all()

class Dotplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("dotplot", *args, **kwargs)
        self.conf.track.load_mat = False

        self.track_a = Track(self.prms.track_a, conf=self.conf)
        self.conf.load_config(self.track_a.conf)

        if self.conf.dotplot.track_b is None:
            self.track_b = None
            self.read_ids = self.track_a.read_ids
        else:
            self.track_b = Track(self.prms.track_b, conf=self.conf)
            self.read_ids = self.track_a.read_ids & self.track_b.read_ids

            if len(self.read_ids) == 0:
                raise RuntimeError("Dotplot tracks must have reads in common")
        
        self.fast5s = Fast5Reader(reads=self.read_ids, conf=self.track_a.conf)

        self.index = self.track_a.index
        self.mm2s = self.track_a.mm2s

    def show_all(self):
        for read_id in self.read_ids:
            self.show(read_id)


    def show(self, read_id, cursor=None):
        if not self._plot(read_id, cursor):
            return False

        print(read_id)

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

    def _plot_aln_step(self, aln, color):
        if getattr(aln, "bands", None) is not None:
            self.ax_dot.fill_between(aln.bands['samp'], aln.bands['ref_st']-1, aln.bands['ref_en'], zorder=1, color='#ccffdd', linewidth=1, edgecolor='black', alpha=0.5)

        self.ax_dot.step(aln.df['start'], aln.df['refmir'], where="post", color=color, zorder=3, linewidth=3)
        #if samp_min is None: samp_min = 0
        #if samp_max is None: samp_max = self.df['sample'].max()
        #i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        model_means = self.track_a.model[aln.df['kmer']]

        pa_diffs = np.abs(aln.df['mean'] - model_means)

        #self.ax_padiff.step(pa_diffs, aln.df['refmir'], color=color, where="post")

    def _plot_signal(self, aln, ax, model):
        samp_min, samp_max = aln.get_samp_bounds()

        samps = np.arange(samp_min, samp_max)
        raw_norm = self.read.get_norm_signal(samp_min, samp_max)

        ymin = np.min(raw_norm[raw_norm>0])
        ymax = np.max(raw_norm[raw_norm>0])
        bases = nt.kmer_base(aln.df['kmer'], 2)

        samp_bases = np.zeros(len(samps), int)
        for i in range(len(aln.df)):
            st = int(aln.df.iloc[i]['start'] - samp_min)
            en = int(st + aln.df.iloc[i]['length'])
            samp_bases[st:en] = bases[i]

        BASE_COLORS = [
            "#80ff80",
            "#8080ff",
            "#ffbd00",
            "#ff8080",
        ]
        for base, color in enumerate(BASE_COLORS):
            ax.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color, interpolate=True)

        ax.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, alpha=0.75, c="#777777")
        ax.step(aln.df['start'], model[aln.df['kmer']], color='white', linewidth=2, where="post")
        ax.vlines(aln.df['start'], ymin, ymax, linewidth=2, color="white")

        evts = (self.read.df['start'] >= samp_min) & (self.read.df['start'] < samp_max) & (self.read.df['norm_sig'] > 0)

        if self.read.has_events:
            ax.step(self.read.df['start'][evts], self.read.df['norm_sig'][evts], where='post', color='black', linewidth=3)
        else:
            ax.scatter(self.read.df['start'][evts], self.read.df['norm_sig'][evts], s=5, alpha=0.75, c="#777777") 

    def set_cursor(self, ref_coord):
        aln,_,_ = self.alns[list(self.focus)[0]]
        refmir = aln.ref_to_refmir(ref_coord)

        i = aln.df['refmir'].searchsorted(refmir)
        samp = aln.df.iloc[i]['start'] + aln.df.iloc[i]['length']/2

        self.cursor = (samp, refmir)

    def _tick_formatter(self, x, pos):
        return self.index.refmir_to_ref(int(x))

    def _load_read(self, read_id):
        self.cursor = None

        fast5_read = self.fast5s[read_id]

        if not read_id in self.mm2s:
            return False

        self.read = ProcRead(fast5_read, conf=self.conf)

        self.aln_bc = BcFast5Aln(self.index, self.read, self.mm2s[read_id], ref_bounds=self.conf.track.ref_bounds)

        self.aln_a = self.track_a.load_aln(read_id, ref_bounds=self.conf.track.ref_bounds)

        if self.track_b is not None and read_id in self.track_b:
            self.aln_b = self.track_b.load_aln(read_id, ref_bounds=self.conf.track.ref_bounds)
        else:
            self.aln_b = None

        self.ref_bounds = self.aln_a.ref_bounds

        return True

    def _plot(self, read_id, cursor=None):
        if not self._load_read(read_id):
            return False

        matplotlib.use("TkAgg")
        plt.style.use(['seaborn'])

        self.fig = plt.figure()
        gspec = self.fig.add_gridspec(
                #ncols=2, nrows=2, 
                ncols=4, nrows=3, 
                width_ratios=[5,1,1,1],
                height_ratios=[1,3,1] 
        )

        self.ax_sig = self.fig.add_subplot(gspec[0,0])
        self.ax_sig2 = self.fig.add_subplot(gspec[2,0])
        self.ax_dot = self.fig.add_subplot(gspec[1,0])
        self.ax_padiff = self.fig.add_subplot(gspec[1,1])
        self.ax_centroid = self.fig.add_subplot(gspec[1,2])
        self.ax_dwell = self.fig.add_subplot(gspec[1,3])

        def sharex(ax1, ax2):
            ax1.get_shared_x_axes().join(ax1,ax2)

        def sharey(ax1, ax2):
            ax1.get_shared_y_axes().join(ax1,ax2)

        sharex(self.ax_dot, self.ax_sig)
        sharex(self.ax_dot, self.ax_sig2)
        sharey(self.ax_dot, self.ax_padiff)
        sharey(self.ax_dot, self.ax_centroid)
        sharey(self.ax_dot, self.ax_dwell)

        fontsize=12
        self.ax_sig.set_ylabel("Current (pA)", fontsize=fontsize)
        self.ax_sig2.set_xlabel("Raw Sample", fontsize=fontsize)
        #self.ax_dot.set_xlabel("Raw Sample", fontsize=fontsize)
        #self.ax_padiff.set_xlabel("Abs. pA Difference", fontsize=fontsize)
        self.ax_padiff.set_xlabel("Sample Jaccard", fontsize=fontsize)
        self.ax_centroid.set_xlabel("Centroid Difference", fontsize=fontsize)
        self.ax_dwell.set_xlabel("Dwell Difference", fontsize=fontsize)

        self.ax_sig.xaxis.set_major_formatter(NullFormatter())
        self.ax_padiff.yaxis.set_major_formatter(NullFormatter())
        self.ax_centroid.yaxis.set_major_formatter(NullFormatter())
        self.ax_dwell.yaxis.set_major_formatter(NullFormatter())

        self.ax_dot.set_ylabel("%s (%s)" % (self.ref_bounds.name, "+" if self.ref_bounds.fwd else "-"), fontsize=12)

        self.ax_dot.yaxis.set_major_formatter(FuncFormatter(self._tick_formatter))
        
        self.ax_sig.set_title(self.read.id)

        if cursor is not None:
            cursor_kw = {
                'color' : 'red', 
                'alpha' : 0.5
            }
            samp, refmir = cursor
            self.ax_dot.axvline(samp, **cursor_kw),
            self.ax_dot.axhline(refmir,  **cursor_kw)

        if self.aln_b is not None:
            self.aln_b.sort_refmir()
            self._plot_signal(self.aln_b, self.ax_sig2, self.track_b.model)
            self._plot_aln_step(self.aln_b, "royalblue")

            compare = method_compare_aln(self.aln_a, self.aln_b)
            self.ax_padiff.step(compare["jaccard"], compare.index, where="post", color="green")
            self.ax_centroid.step(compare["centroid_diff"], compare.index, where="post", color="green")
            self.ax_dwell.step(compare["dwell_diff"], compare.index, where="post", color="green")

        self.aln_a.sort_refmir()
        self._plot_signal(self.aln_a, self.ax_sig, self.track_a.model)
        self._plot_aln_step(self.aln_a, "purple")

        self.ax_dot.scatter(self.aln_bc.df['sample'], self.aln_bc.df['refmir'], color='orange', zorder=2, s=20)

        self.fig.tight_layout()

        return True


    def save(self, out_prefix, fmt):
        self._plot()

        suffix = self.read.id + "." + fmt
        if os.path.isdir(out_prefix):
            out_fname = os.path.join(out_prefix, suffix)
        else:
            out_fname = out_prefix + suffix

        self.fig.set_size_inches(8, 6)
        self.fig.savefig(out_fname, dpi=200)
        plt.close()
