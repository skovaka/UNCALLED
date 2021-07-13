import time
import numpy as np
import os 

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter
import types

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
    ("layers", ["current", "length"], list, "Layers to plot in side panels"),
    ("style", {
        "aln_kws" : [
            {"color" : "purple", "alpha" : 1},
            {"color" : "royalblue", "alpha" : 0.5}],
        "bc_color": "orange",
        "base_colors" : ["#80ff80", "#8080ff", "#ffbd00", "#ff8080"], #A,C,G,T
    }, dict, "Plotting style options")
)
#Config._EXTRA_GROUPS["dotplot"] = DotplotParams #TODO put in ParamGroup con

def comma_split(s):
    return s.split(",")

from ..config import Opt
OPTS = [
    Opt("track_a", "dotplot"),
    Opt("track_b", "dotplot", nargs="?"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "track", type=RefCoord),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", type=comma_split),
    Opt(("-C", "--max-chunks"), "read_buffer"),
]

def main(conf):
    """Plot dotplots of alignments from tracks produced by `align` or `convert`"""
    Dotplot(conf=conf).show_all()

class Dotplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("dotplot", *args, **kwargs)
        self.conf.track.load_mat = False
        self.conf.track.layers = self.conf.dotplot.layers

        self.tracks = [Track(self.prms.track_a, conf=self.conf)]
        #self.colors = [self.prms.styles["color_a"]]

        self.conf.load_config(self.tracks[0].conf)

        if self.conf.dotplot.track_b is not None:
            self.tracks.append(Track(self.prms.track_b, conf=self.conf))
            #self.colors.append(self.prms.styles["color_b"])

            self.read_ids = self.tracks[0].read_ids & self.tracks[1].read_ids

            if len(self.read_ids) == 0:
                raise RuntimeError("Dotplot tracks must have reads in common")
        else:
            self.read_ids = self.tracks[0].read_ids
        
        self.fast5s = Fast5Reader(reads=self.read_ids, conf=self.tracks[0].conf)

        self.index = self.tracks[0].index
        self.mm2s = self.tracks[0].mm2s

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

    def _plot_signal(self, aln, ax, model):
        samp_min, samp_max = aln.get_samp_bounds()

        samps = np.arange(samp_min, samp_max)
        raw_norm = self.read.get_norm_signal(samp_min, samp_max)

        ymin = np.min(raw_norm[raw_norm>0])
        ymax = np.max(raw_norm[raw_norm>0])
        bases = nt.kmer_base(aln.aln['kmer'], 2)

        samp_bases = np.zeros(len(samps), int)
        for i in range(len(aln.aln)):
            st = int(aln.aln.iloc[i]['start'] - samp_min)
            en = int(st + aln.aln.iloc[i]['length'])
            samp_bases[st:en] = bases[i]

        for base, color in enumerate(self.prms.style["base_colors"]):
            ax.fill_between(samps, ymin, ymax, where=samp_bases==base, color=color, interpolate=True)

        ax.scatter(samps[raw_norm > 0], raw_norm[raw_norm > 0], s=5, alpha=0.75, c="#777777")
        ax.step(aln.aln['start'], model[aln.aln['kmer']], color='white', linewidth=2, where="post")
        ax.vlines(aln.aln['start'], ymin, ymax, linewidth=2, color="white")

        evts = (self.read.df['start'] >= samp_min) & (self.read.df['start'] < samp_max) & (self.read.df['norm_sig'] > 0)

        if self.read.has_events:
            ax.step(self.read.df['start'][evts], self.read.df['norm_sig'][evts], where='post', color='black', linewidth=3)
        else:
            ax.scatter(self.read.df['start'][evts], self.read.df['norm_sig'][evts], s=5, alpha=0.75, c="#777777") 

    def _plot_aln(self, i):
        aln = self.alns[i]
        if not "refmir" in aln.aln:
            aln.calc_refmir()
        aln.sort_refmir()

        if getattr(aln, "bands", None) is not None:
            self.ax_dot.fill_between(aln.bands['samp'], aln.bands['ref_st']-1, aln.bands['ref_en'], zorder=1, color='#ccffdd', linewidth=1, edgecolor='black', alpha=0.5)

        self.ax_dot.step(aln.aln['start'], aln.aln['refmir'], where="post", linewidth=3,
            **self.prms.style["aln_kws"][i]
        )

        #if "bcerr" in aln.dfs:
        #    print(aln.ref_to_samp(aln.bcerr[aln.bcerr["type"] == "INS"].index))

        self._plot_signal(aln, self.sig_axs[i], self.tracks[i].model)

        #if samp_min is None: samp_min = 0
        #if samp_max is None: samp_max = self.df['sample'].max()
        #i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)
        #model_means = self.track_a.model[aln.aln['kmer']]
        #pa_diffs = np.abs(aln.aln['mean'] - model_means)
        #self.ax_padiff.step(pa_diffs, aln.aln['refmir'], color=color, where="post")

    def set_cursor(self, ref_coord):
        aln,_,_ = self.alns[list(self.focus)[0]]
        refmir = aln.ref_to_refmir(ref_coord)

        i = aln.aln['refmir'].searchsorted(refmir)
        samp = aln.aln.iloc[i]['start'] + aln.aln.iloc[i]['length']/2

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

        self.alns = [
            track.load_read(read_id)#, ref_bounds=self.conf.track.ref_bounds)
            for track in self.tracks
        ]

        #TODO improve this
        self.ref_bounds = self.alns[0].ref_bounds

        return True
    
    #def _load_layer(self, layer):

    
    @property
    def track_count(self):
        return len(self.tracks)

    def _plot(self, read_id, cursor=None):
        if not self._load_read(read_id):
            return False

        matplotlib.use("TkAgg")
        plt.style.use(['seaborn'])

        widths=[5]
        for _ in self.prms.layers:
            widths.append(1)

        heights=[1,3] 
        if self.track_count == 2:
            heights.append(1)

        self.fig = plt.figure()
        gspec = self.fig.add_gridspec(
                ncols=len(widths), nrows=len(heights), 
                width_ratios=widths,
                height_ratios=heights 
        )

        self.sig_axs = [self.fig.add_subplot(gspec[i*2,0]) 
                        for i in range(len(self.tracks))]

        self.ax_dot = self.fig.add_subplot(gspec[1,0])

        self.layer_axs = [self.fig.add_subplot(gspec[1,l+1])
                          for l in range(len(self.prms.layers))]

        fontsize=12

        for ax in self.sig_axs:
            ax.get_shared_x_axes().join(self.ax_dot, ax)
            ax.set_ylabel("Current (pA)", fontsize=fontsize)

        self.sig_axs[0].set_title(self.read.id)
        self.sig_axs[-1].set_xlabel("Raw Sample", fontsize=fontsize)

        for i,ax in enumerate(self.layer_axs):
            ax.get_shared_y_axes().join(self.ax_dot, ax)
            ax.yaxis.set_major_formatter(NullFormatter())
            ax.set_xlabel(self.prms.layers[i], fontsize=fontsize)

        self.ax_dot.set_ylabel("%s (%s)" % (self.ref_bounds.name, "+" if self.ref_bounds.fwd else "-"), fontsize=12)

        self.ax_dot.yaxis.set_major_formatter(FuncFormatter(self._tick_formatter))
        
        if cursor is not None:
            cursor_kw = {
                'color' : 'red', 
                'alpha' : 0.5
            }
            samp, refmir = cursor
            self.ax_dot.axvline(samp, **cursor_kw),
            self.ax_dot.axhline(refmir,  **cursor_kw)

        for i in range(self.track_count):
            self._plot_aln(i)

        for ax,layer in zip(self.layer_axs, self.prms.layers):
            for i,aln in enumerate(self.alns):
                ax.step(
                    aln.aln[layer], aln.aln["refmir"], 
                    where="post",
                    **self.prms.style["aln_kws"][i]
                )


        if not self.aln_bc.empty:
            self.ax_dot.scatter(self.aln_bc.aln['sample'], self.aln_bc.aln['refmir'], color='orange', zorder=2, s=20)

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
