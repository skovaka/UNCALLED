import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time
import sys

from .. import config
from ..dtw.layers import LAYER_META, parse_layer, parse_layers
from ..index import str_to_coord
from ..dtw.tracks import Tracks, REFSTAT_LABELS, COMPARE_REFSTATS
from ..argparse import Opt, comma_split

class Refplot:

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("refplot", *args, **kwargs)

        if self.prms.tracks is None:
            self.tracks = Tracks(conf=self.conf)
        else:
            self.tracks = self.prms.tracks
            self.tracks.conf.load_config(self.conf)
            self.conf = self.tracks.conf

        if self.tracks.refstats is None:
            self.tracks.calc_refstats()

        names = [t.name for t in self.tracks.alns]

        if self.prms.kmer_coord is None:
            self.refs = self.tracks.coords.refs
            mid_plot = len(self.refs)//2
        else:
            st = en = self.prms.kmer_coord 
            st -= self.tracks.index.trim[0]
            en += self.tracks.index.trim[1] + 1
            self.refs = pd.RangeIndex(st, en)
            mid_plot = self.tracks.index.trim[0]

        self.layer, = parse_layer(self.prms.layer)

        subplot_titles = list()
        if self.prms.kmer_coord is not None:
            for i in range(len(self.refs)):
                shift = i-mid_plot
                if shift < 0:
                    subplot_titles.append(f"{shift}")
                elif shift == 0:
                    subplot_titles.append("Selection")
                else:
                    subplot_titles.append(f"+{shift}")

        self.fig = make_subplots(
            rows=1, cols=len(self.refs), 
            horizontal_spacing=0.025,
            subplot_titles=subplot_titles,
            shared_xaxes=True, 
            shared_yaxes=True)
            #vertical_spacing=0.125/n_rows)

        self.fig.update_layout(
            dragmode="pan", 
            height=300,
            margin=dict(l=65,r=25,b=60,t=55),
            showlegend=False,
        )

        self.fig.update_xaxes(title=self.tracks.coords.ref_name, row=1, col=mid_plot+1)
            #title=self.tracks.coords.ref_name, title_x=0.5, title_y=0.05,

        #if self.prms.share_reads:
        self.fig.update_yaxes(
            title={"text" : LAYER_META.loc[self.layer,"label"], "standoff" : 8},
            row=1,col=1)
        self.fig.update_xaxes(matches="x", tickvals=[0], tickmode="array")


        legend = set()

        for r,ref in enumerate(self.refs):
            for i,track in enumerate(self.tracks.alns):
                vals = self.tracks.layers.loc[(track.id,slice(None),ref),self.layer].dropna()
                counts, bins = np.histogram(vals, bins=25)
                dens = counts / counts.sum()
                if len(self.tracks.alns) == 2 and i == 1:
                    dens = -dens
                width = bins[1]-bins[0]

                showlegend = track.name not in legend
                if showlegend:
                    legend.add(track.name)

                self.fig.add_trace(
                    go.Bar(
                        y=bins[1:]-width/2, x=dens, 
                        orientation="h",
                        width=width,
                        legendgroup=track.name,
                        showlegend=showlegend,
                        name=track.desc,
                        marker={
                            "line_width" : 0,
                            "color" : self.conf.vis.track_colors[i],
                    }), row=1, col=r+1
                )
                self.fig.update_xaxes(ticktext=[str(ref)], row=1, col=r+1)
                

    def show(self):
        fig_conf = {
            "toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None, "scale" : 2},
            "scrollZoom" : True, 
            "displayModeBar" : True}
        self.fig.show(config=fig_conf)

def refplot(conf):
    """Plot per-reference distributions"""
    conf.tracks.layers.append(conf.refplot.layer)
    Refplot(conf=conf).show()
