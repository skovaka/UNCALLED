import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time
import sys

from .. import config
from ..fast5 import parse_read_ids
from ..dtw.aln_track import LAYERS, parse_layer, parse_layers
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

        names = [t.name for t in self.tracks]

        if self.prms.kmer_coord is None:
            self.refs = self.tracks.coords.refs
        else:
            st = en = self.prms.kmer_coord 
            st -= self.tracks.index.trim[0]
            en += self.tracks.index.trim[1] + 1
            self.refs = pd.RangeIndex(st, en)

        self.layer, = parse_layer(self.prms.layer)

        self.fig = make_subplots(
            rows=1, cols=len(self.refs), 
            shared_xaxes=True, 
            shared_yaxes=True)
            #vertical_spacing=0.125/n_rows)

        self.fig.update_layout(
            title=self.tracks.coords.ref_name, title_x=0.5, title_y=0.05,
            height=300,
            margin=dict(l=75,r=25,b=75,t=25),
            showlegend=False,
        )

        #if self.prms.share_reads:
        self.fig.update_yaxes(title=LAYERS[self.layer[0]][self.layer[1]].label,row=1,col=1)
        self.fig.update_xaxes(matches="x", showticklabels=False)

        legend = set()

        for r,ref in enumerate(self.refs):
            for i,track in enumerate(self.tracks.alns):
                vals = track.layers.loc[ref,self.layer]
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
                self.fig.update_xaxes(title=str(ref), row=1, col=r+1)

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
