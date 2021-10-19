import sys, os
import numpy as np
import argparse
from collections import defaultdict
import re
import time
import types

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter
from matplotlib.colors import Normalize, TwoSlopeNorm, ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib import widgets
import matplotlib

import scipy.stats
import pandas as pd

from ... import config, nt

#from ..config import Config, ParamGroup, Opt
from ...sigproc import ProcRead

from ...index import BWA_OPTS, str_to_coord
from ...fast5 import Fast5Reader
from ..dotplot import Dotplot
from ...dtw.track import AlnTrack, LAYERS
from ...dtw.tracks import Tracks
from ...argparse import Opt, comma_split

CMAP = "viridis"
#CMAP = "plasma"

class BrowserParams(config.ParamGroup):
    _name = "browser"
BrowserParams._def_params(
    ("track_a", None, None, "Path to directory where alignments are stored"),
    ("track_b", None, None, "Path to directory where alignments are stored"),
    ("full_overlap", None, bool, "If true will only include reads which fully cover reference bounds"),
    ("style", {
        "rc" : {
            "figure.figsize" : (16, 10)
        },
    }, dict, "Plotting style options")
)

OPTS = (
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt("input", "tracks", nargs="+"),
    #Opt("track_a", "browser", type=str),
    #Opt("track_b", "browser", type=str, nargs="?"),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-L", "--layers"), "tracks", type=comma_split),
)

def main(conf):
    """Display interactive signal-based genome browser"""

    matplotlib.use("TkAgg")
    browser = Browser(conf=conf)
    browser.show()

class Browser:

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("browser", *args, **kwargs)

        matplotlib.rcdefaults()
        plt.style.use(['seaborn'])
        matplotlib.rcParams.update(self.prms.style["rc"])

        self.io = Tracks(conf=self.conf)
        self.conf = self.io.conf

        self.tracks = self.io.load_refs(load_mat=True) #list()

        self.single_track = len(self.tracks) == 1

        ref_bounds = self.io.prms.ref_bounds

        self.ref_name = ref_bounds.name
        self.ref_start = ref_bounds.start
        self.ref_end = ref_bounds.end
        self.width = self.ref_end - self.ref_start


        if not self.single_track:
            self.tracks = self.tracks[:2]
            self.ks_stats = self.tracks[0].calc_ks(self.tracks[1])

        #self.fast5s = #Fast5Reader(conf=self.conf)
        self.fast5s = self.io.get_fast5_reader()#Fast5Reader(reads=self.read_ids, index=fast5_index, conf=self.tracks[0].conf)

        self.cursor = None
        self.pileup = False

        self.LAYER_IDS = dict()
        for i,layer in enumerate(self.conf.tracks.layers):
            self.LAYER_IDS[LAYERS[layer].label] = layer

        self.active_layer = self.conf.tracks.layers[0] #TODO parameter

        self.info_labels = ["Reference Coord", "Read ID"]
        for layer in self.conf.tracks.layers:
            self.info_labels.append(LAYERS[layer].label)

        self.axs = types.SimpleNamespace()

        self.shown = False
    
    @property
    def ref_coords(self):
        #TODO store and/or make sure they match
        return "%s:%d-%d" % (self.ref_name, self.ref_start, self.ref_end)
    
    @property
    def cursor_coords(self):
        return self.cursor[0].get_xdata(), self.cursor[1].get_ydata()

    def plot_dotplot(self, event):
        tr,rf,rd = self.cursor
        track = self.tracks[tr]

        read = track.alignments.iloc[rd]

        #fast5_read = self.fast5s[read.id]
        #proc_read = ProcRead(fast5_read, conf=self.conf)
        #aln = track.load_aln(read.id)
        #bcaln = BcFast5Aln(aln.index, proc_read, track.mm2s[read.id])

        dpl = Dotplot(track, conf=self.conf)
        dpl.show(read.id)
        #if not bcaln.empty:
        #    dpl.add_aln(bcaln, False)
        #dpl.add_aln(aln, True)
        #dpl.set_cursor(self.ref_start + rf)
        #dpl.show()


        #fig = plt.figure()
        #ax = fig.subplots(1,1)
        #dtw.plot_dotplot(ax)
        #fig.show()

    def toggle_pileup(self, event=None):
        self.pileup = not self.pileup
        for track in self.tracks:
            track.img.set_data(self.get_img(track))
        self.fig.canvas.draw()
        #if self.pileup:

    def sort_position(self, event=None):
        tr,rf,rd = self.cursor
        self.pileup = False
        for track in self.tracks:
            track.sort_coord(self.active_layer, self.ref_start + rf)
            track.img.set_data(self.get_img(track))
            self.del_cursor()

        self.fig.canvas.draw()

    def down_click(self, event):
        self.dclick_x = event.x
        self.dclick_y = event.y

    def track_ax_index(self, ax):
        for i,trk in enumerate(self.tracks):
            if trk.ax == ax:
                return i
        return -1


    def up_click(self,event):
        tr = self.track_ax_index(event.inaxes)
        if (tr < 0 or 
            event.button != 1 or
            abs(self.dclick_x - event.x) > 0.1 or 
            abs(self.dclick_y - event.y) > 0.1):  #TODO constant

            return 
        
        track = self.tracks[tr]
        rf, rd = int(event.xdata+0.5), int(event.ydata+0.5)

        self.set_cursor(tr, rf, rd)

        #if self.single_track:
        self.set_info(track, rf, rd)
        #self.plot_hists(track, rf)

        self.fig.tight_layout()
        self.fig.canvas.draw()

    def del_cursor(self):
        self.cursor = None
        to_hide = [self.axs.info, self.axs.pa_hist, self.axs.dwell_hist, self.axs.dot_btn, self.axs.sort_btn]
        for track in self.tracks:
            to_hide += track.cursor

        for ax in to_hide:
            ax.set_visible(False)

    def set_cursor(self, tr, rf, rd):
        self.cursor = (tr, rf, rd)

        for i,trk in enumerate(self.tracks):
            trk.cursor[0].set_xdata(rf)
            trk.cursor[0].set_visible(True)

            if i == tr:
                trk.cursor[1].set_ydata(rd)
                trk.cursor[1].set_visible(True)
            else:
                trk.cursor[1].set_visible(False)

    def set_info_cell(self, row, val):
        self.info_table[row, 1].get_text().set_text(val)

    def ref_coord(self, track, i, fwd=None):
        coord = "%s:%d" % (
                self.ref_name, 
                self.ref_start + i
        )
        if fwd is not None:
            coord += " (%s)" % ("+" if fwd else "-")
        return coord

    def set_info(self, trk, rf, rd):

        aln_id = trk.alignments.index[rd]

        read = trk.alignments.iloc[rd]
        #aln_id = trk.aln

        self.set_info_cell(0,  self.ref_coord(trk, rf, read['fwd']))
        self.set_info_cell(1, read['read_id'])

        mref = trk.coords.ref_to_mref([rf], True)[0]

        for i,layer in enumerate(self.conf.tracks.layers):
            val = trk.layers["dtw",layer].loc[mref,aln_id]
            self.set_info_cell(i+2, val)
        #self.set_info_cell(self.INFO_KMER, nt.kmer_to_str(int(kmer)))
        #self.set_info_cell(self.INFO_PA,   "%.4f" % pa)
        #self.set_info_cell(self.INFO_DWELL,"%.4f" % dwell)
        #self.set_info_cell(self.INFO_DIFF, "%.4f" % diff)

        self.info_table.auto_set_column_width(1)

        self.axs.info.set_visible(True)
        self.axs.pa_hist.set_visible(True)
        self.axs.dwell_hist.set_visible(True)
        self.axs.dot_btn.set_visible(True)
        self.axs.sort_btn.set_visible(True)

    
    def plot_hists(self, track, rf):
        self.axs.pa_hist.cla()
        self.axs.dwell_hist.cla()

        fwds = track.alignments['fwd']

        pa    = track[self.PA_LAYER,   :,rf]
        pa_fwd = pa[fwds]
        pa_rev = pa[~fwds]

        dwell = track[self.DWELL_LAYER,:,rf]
        dwell_fwd = dwell[fwds]
        dwell_rev = dwell[~fwds]

        kw = {
            "density" : True,
            "alpha"   : 0.75
        }
        
        fwd_kw = {"color" : "royalblue"}
        rev_kw = {"color" : "crimson"}
        fwd_kw.update(kw)
        rev_kw.update(kw)

        pa_min = np.min(pa)
        pa_max = np.max(pa)
        pa_span = pa_max-pa_min

        pa_bins = int(pa_span) #TODO constant
        dwell_bins = int((np.max(dwell)-np.min(dwell))/2) #TODO constant


        if track.has_fwd:
            self.axs.pa_hist.hist(pa_fwd, bins=pa_bins, **fwd_kw)
            self.axs.dwell_hist.hist(dwell_fwd, bins=dwell_bins, **fwd_kw)

            kmer = int(scipy.stats.mode(track[self.KMER_LAYER,track.alignments['fwd'],rf]).mode[0])
            exp_mean = track.model.means[kmer]
            exp_stdv = track.model.stdvs[kmer]

            slop = pa_span*0.05
            xs = np.linspace(pa_min-slop, pa_max+slop, 100)
            ys = np.exp(track.model.norm_pdf(xs, kmer))
            #ys = scipy.stats.norm.pdf(xs, exp_mean, exp_stdv)
            self.axs.pa_hist.plot(xs, ys, color="blue", linewidth=2)

        if track.has_rev:
            self.axs.pa_hist.hist(pa_rev, bins=pa_bins, **rev_kw)
            self.axs.dwell_hist.hist(dwell_rev, bins=dwell_bins, **rev_kw)

            kmer = int(scipy.stats.mode(track[self.KMER_LAYER,~fwds,rf]).mode[0])
            exp_mean = track.model.means[kmer]
            exp_stdv = track.model.stdvs[kmer]
            slop = pa_span*0.05
            xs = np.linspace(pa_min-slop, pa_max+slop, 100)
            ys = np.exp(track.model.norm_pdf(xs, kmer))
            #ys = scipy.stats.norm.pdf(xs, exp_mean, exp_stdv)
            self.axs.pa_hist.plot(xs, ys, color="red", linewidth=2)

        self.axs.pa_hist.set_xlabel(self.LAYERS[self.PA_LAYER].label)
        self.axs.dwell_hist.set_xlabel(self.LAYERS[self.DWELL_LAYER].label)

        self.axs.pa_hist.set_yticks([])
        self.axs.dwell_hist.set_yticks([])

        #hist_fig.tight_layout()
        #hist_fig.show()

    def set_layer(self, layer):
        if layer in self.LAYER_IDS:
            layer = self.LAYER_IDS[layer]

        if layer == self.active_layer: return

        self.active_layer = layer

        if not self.shown: return

        for track in self.tracks:

            track.img.set_data(self.get_img(track))

            cmap, norm = self._layer_cmap(self.active_layer)

            #TODO same norm between images
            track.img.set_norm(norm)
            track.img.set_cmap(cmap)

        #self.cbar.update_normal(ScalarMappable(norm=track.norms[layer], cmap=CMAP))
        #self.cbar.update_normal(ScalarMappable(norm=norm, cmap=cmap))
        self.cbar.update_normal(ScalarMappable(norm=norm, cmap=cmap))
        self.axs.cbar.set_ylabel(LAYERS[layer].label) 
        self.axs.cbar.yaxis.set_label_position('left') 

        self._update_sumstat()

        self.fig.canvas.draw()

    def _update_sumstat(self):
        if self.single_track:
            track = self.tracks[0]
            if track.has_fwd:
                self.layer_fwd_plot.set_ydata(self.layer_means(track, self.active_layer, True))

            if track.has_rev:
                self.layer_rev_plot.set_ydata(self.layer_means(track, self.active_layer, False))

            self.axs.sumstat.relim()
            self.axs.sumstat.autoscale_view()

        #else:
        #    stats = self.calc_ks(self.active_layer)
        #    self.sumstat_img.set_data(stats)


    def ref_tick_fmt(self, x, pos=None):
        return int(np.round(self.ref_start + x))
    

    #TODO put inside track 
    def layer_means(self, track, layer, fwd):
        #return np.mean(track.mat[layer], axis=0)
        x = track.mat["dtw"][layer].mean(axis=0)#.reindex(track.coords.refs)
        return x

    #TODO refactor into _init_fig(), show(), and save()
    def show(self):
        self.fig = plt.figure(
            tight_layout={"pad" : 2, "h_pad" : 0, "w_pad" : 1}
        )

        gs_main = self.fig.add_gridspec(
                ncols=2, nrows=1, 
                width_ratios=[5,2] 
        )

        if self.single_track:
            self._init_single_track(gs_main)
        else:
            self._init_double_track(gs_main)

        for track in self.tracks:
            track.ax.get_shared_x_axes().join(track.ax, self.axs.sumstat)
        
        self.fig.canvas.mpl_connect("close_event", self.on_close)
        self.fig.canvas.mpl_connect('button_press_event', self.down_click) 
        self.fig.canvas.mpl_connect('button_release_event', self.up_click) 

        self.shown = True
        #mng = plt.get_current_fig_manager()
        #mng.window.attributes('-zoomed', True)

        plt.show()

    def _init_single_track(self, gspec):
        gs_left = gspec[0].subgridspec(2, 1, height_ratios=[4,1], hspace=0)

        self._init_track_ax(self.tracks[0], gs_left[0])
        self._init_sumstat_mean(gs_left[1], self.tracks[0])

        gs_right = gspec[1].subgridspec(5, 1, height_ratios=[4,5,4,4,1], hspace=0.35)
        self._init_layer_picker(gs_right[0])
        self._init_info(gs_right[1])
        self._init_hists(gs_right[2], gs_right[3])
        self._init_btns(gs_right[4])

    def _init_double_track(self, gspec):
        gs_left = gspec[0].subgridspec(3, 1, height_ratios=[2,2,1], hspace=0.05)
        self._init_track_ax(self.tracks[0], gs_left[0])
        self._init_track_ax(self.tracks[1], gs_left[1], False)
        self._init_sumstat_ks(gs_left[2])


        #gs_right = gspec[1].subgridspec(2, 1, height_ratios=[4,14], hspace=0.35)
        gs_right = gspec[1].subgridspec(5, 1, height_ratios=[4,5,4,4,1], hspace=0.35)
        self._init_layer_picker(gs_right[0])
        self._init_info(gs_right[1])
        self._init_hists(gs_right[2], gs_right[3])
        self._init_btns(gs_right[4])
    
    @staticmethod
    def _noticks(ax, x=True, y=True):
        if x:
            ax.xaxis.set_major_formatter(NullFormatter())
        if y:
            ax.yaxis.set_major_formatter(NullFormatter())

    @staticmethod
    def _nospines(ax):
        for name,spine in ax.spines.items():
            spine.set_visible(False)

    def _init_layer_picker(self, gspec):
        subspec = gspec.subgridspec(1, 2, width_ratios=[1,6], wspace=0.1)

        self.axs.cbar = self.fig.add_subplot(subspec[0])
        cmap, norm = self._layer_cmap(self.active_layer)
        self.cbar = self.fig.colorbar(
            ScalarMappable(norm=norm, cmap=cmap),
            cax=self.axs.cbar, extend="both", orientation="vertical"
        )
        self.cbar.outline.set_visible(False)
        self.axs.cbar.yaxis.set_label_position('left') 
        self.axs.cbar.set_ylabel(LAYERS[self.active_layer].label) 

        self.axs.opt = self.fig.add_subplot(subspec[1])

        self.axs.opt.patch.set_alpha(0.0)
        self._nospines(self.axs.opt)
        self._noticks(self.axs.opt)

        self.layer_radio = widgets.RadioButtons(
            self.axs.opt, [LAYERS[layer].label for layer in self.conf.tracks.layers],
            0, #TODO hacky way to deal with hidden layers
            activecolor = 'red'
        )
        self.layer_radio.on_clicked(self.set_layer)

    #LAYER_CMAPS = {
    #    "model_diff" : lambda track: ("bwr", TwoSlopeNorm(vcenter=0))
    #}

    def _layer_cmap(self, layer):
        if layer == "bcerr":
            return (ListedColormap(["#80ff80", "#8080ff", "#ffbd00", "#ff8080", "purple", "yellow"]), None)

        mins = list()
        maxs = list()
        for track in self.tracks:
            std = track.layers["dtw"][layer].std()
            med = track.layers["dtw"][layer].median()
            mins.append(min(track.layers["dtw"][layer].min(), med - std*3))
            maxs.append(min(track.layers["dtw"][layer].max(), med + std*3))

        vmin = min(mins)
        vmax = max(maxs)

        if layer == "model_diff":
            return ("bwr", TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax))
        return ("viridis", Normalize(vmin, vmax))

    def _init_track_ax(self, track, gspec, ref_ticks=True):
        track.ax = self.fig.add_subplot(gspec) 

        if ref_ticks:
            track.ax.xaxis.set_label_position('top') 
            track.ax.xaxis.tick_top()
            track.ax.xaxis.set_major_formatter(FuncFormatter(self.ref_tick_fmt))

        track.ax.grid()
        self._noticks(track.ax, x=not ref_ticks)

        track.ax.set_xlabel(self.ref_coords, fontsize=15)
        track.ax.set_ylabel("%s (%d reads)" % (track.name, track.height), fontsize=15)

        cmap, norm = self._layer_cmap(self.active_layer)
        track.img = track.ax.imshow(
            self.get_img(track), 
            #alpha=np.array((~track.track.mask[0]).astype(float)), #TODO use mat_alpha()
            aspect='auto', 
            cmap=cmap, 
            norm=norm,
            interpolation='none'#'antialiased' if self.zoomed_out() else 'none'
        )

        self._init_cursor(track)

    def get_img(self, track):
        if self.pileup:
            return track.get_pileup(self.active_layer)
        return track.mat["dtw"][self.active_layer]

    def _init_cursor(self, track):
        cursor_kw = {
            'color' : 'red', 
            'alpha' : 0.75,
            'visible' : self.cursor != None
        }

        track.cursor = (
            track.ax.axvline(0,0,track.height, **cursor_kw),
            track.ax.axhline(0,0,track.width,  **cursor_kw)
        )

        if self.cursor != None:
            self.set_cursor(*self.cursor)
            

    def _init_sumstat(self, gspec):
        self.axs.sumstat = self.fig.add_subplot(gspec)
        self.axs.sumstat.xaxis.set_major_formatter(FuncFormatter(self.ref_tick_fmt))

    def _init_sumstat_ks(self, gspec):
        self._init_sumstat(gspec)
        ax = self.axs.sumstat
        ax.grid()

        self.sumstat_img = ax.imshow(
            self.ks_stats, 
            aspect="auto", 
            cmap="inferno", 
            interpolation='none'
        )

        ax.set_yticks(np.arange(len(AlnTrack.CMP_LAYERS)))
        ax.set_yticklabels([LAYERS[l].label.split()[0] for l in AlnTrack.CMP_LAYERS])
            
    def _init_sumstat_mean(self, gspec, track):
        self._init_sumstat(gspec)
        self.axs.sumstat.set_ylabel("Mean Value", fontsize=15) #TODO use rcparams

        if track.has_fwd:
            self.layer_fwd_plot, = self.axs.sumstat.plot(np.arange(track.width), self.layer_means(track, self.active_layer, True), color='royalblue')

        if track.has_rev:
            self.layer_rev_plot, = self.axs.sumstat.plot(np.arange(track.width), self.layer_means(track, self.active_layer, False), color='crimson')

    def _init_info(self, gspec):
        self.axs.info  = self.fig.add_subplot(gspec)
        self._nospines(self.axs.info)
        self.axs.info.set_xticks([])
        self.axs.info.set_yticks([])
        self.axs.info.set_facecolor('white')

        y = 2*len(self.info_labels)
        self.axs.info.set_ylim(0, y+0.5)
        self.axs.info.set_xlim(-0.1, 10)

        self.info_table = self.axs.info.table(
            cellText = [[label, ""] for label in self.info_labels],
            cellColours = [["lightcoral", "white"] for _ in self.info_labels],
            loc='upper center', cellLoc='left'
        )
        self.info_table.auto_set_font_size(True)
        
        for i in range(len(self.info_labels)):
            self.info_table[i,1].PAD=0.05
            self.info_table[i,0].set_text_props(weight='bold')
            self.info_table[i,1].set_text_props(family='monospace')

        self.info_table.auto_set_column_width(0)
        self.axs.info.set_visible(False)

    def _init_hists(self, gs_pa, gs_dwell):
        self.axs.pa_hist    = self.fig.add_subplot(gs_pa)
        self.axs.dwell_hist = self.fig.add_subplot(gs_dwell)
        self.axs.pa_hist.set_visible(False)
        self.axs.dwell_hist.set_visible(False)

    def _init_btns(self, gspec):
        subspec = gspec.subgridspec(1, 3, wspace=0.1)

        self.axs.pileup_btn = self.fig.add_subplot(subspec[0])
        self.pileup_btn = widgets.Button(self.axs.pileup_btn, "Pileup")
        self.pileup_btn.on_clicked(self.toggle_pileup)

        self.axs.sort_btn = self.fig.add_subplot(subspec[1])
        self.sort_btn = widgets.Button(self.axs.sort_btn, "Sort")
        self.sort_btn.on_clicked(self.sort_position)
        self.axs.sort_btn.set_visible(False)

        self.axs.dot_btn = self.fig.add_subplot(subspec[2])
        self.dot_btn = widgets.Button(self.axs.dot_btn, "Dotplot")
        self.dot_btn.on_clicked(self.plot_dotplot)
        self.axs.dot_btn.set_visible(False)

    def on_close(self, event):
        self.io.close()
        #for t in self.tracks:
        #    t.close()
        sys.exit(0)


