import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time

from ... import config
from ...dtw.track import LAYERS
from ...index import str_to_coord
from ...dtw.track_io import Tracks
from ...argparse import Opt, comma_split

class TrackplotParams(config.ParamGroup):
    _name = "trackplot"
TrackplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("layer", "current", str, "Layer to plot"),
    ("track_colors", ["#AA0DFE", "#1CA71C", "#6A76FC"], list, ""),
    ("select_ref", None, str, "Reference Selection"),
    ("select_read", None, str, "Read Selection"),
    ("outfile", None, str, "Output file"),
    #("track_colors", ["#AA0DFE", "#1CA71C", "#6A76FC"], list, ""),
)

LAYER_COLORS = {
    "model_diff" : {"colorscale" : "RdBu", "cmid" : 0, "cmax" : 20, "cmin" : -20, "reversescale":True},
    "current" : {"colorscale" : "viridis"},
    "dwell" : {"colorscale" : "viridis", "cmin" : 0, "cmax" : 25},
}

class Trackplot:

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("trackplot", *args, **kwargs)

        if self.prms.tracks is None:
            self.tracks = Tracks(conf=self.conf)
            self.tracks.load_refs(load_mat=True)
        else:
            self.tracks = self.prms.tracks

        names = [t.name for t in self.tracks]

        if self.tracks.refstats is None:
            refstat_tracks = pd.Index([])
            layer_stats = pd.Index([])
        else:
            refstat_tracks = self.tracks.refstats.columns.get_level_values("track").unique()
            cmp_stats = refstat_tracks.difference(names)
            if len(refstat_tracks.intersection(names)) == len(names):
                layer_stats = self.tracks.refstats[names].columns.get_level_values("stat").unique()
            else:
                layer_stats = pd.Index([])

        row_heights = [4]*len(self.tracks) + [1]*(len(layer_stats)+len(cmp_stats))
        n_rows = len(row_heights)

        t0 = time.time()
        self.fig = make_subplots(
            rows=n_rows, cols=1, 
            #subplot_titles=names, 
            row_heights=row_heights,
            shared_xaxes=True, 
            x_title=self.tracks.coords.ref_name,
            #y_title="Reads",
            vertical_spacing=0.125/n_rows)

        layer_label = LAYERS["dtw"][self.prms.layer].label

        t0 = time.time()
        for i,track in enumerate(self.tracks.all):
            mat = track.mat["dtw"][self.prms.layer]
            hover = "<br>".join([
                track.coords.ref_name + ":%{x}",
                #"Read: %{y}",
                layer_label + ": %{z}"])

            self.fig.add_trace(go.Heatmap(
                name=track.name,
                x=mat.columns,
                #y=track.alignments["read_id"],
                z=mat,
                zsmooth=False,
                hovertemplate=hover,
                coloraxis="coloraxis",
            ), row=i+1, col=1)

            self.fig.update_yaxes(
                title_text="<b>"+track.name+"</b>", 
                title_font_color=self.prms.track_colors[i], 
                row=i+1, col=1)

            if self.prms.select_read is not None:
                ys = np.where((self.prms.select_read == track.alignments["read_id"]).to_numpy())[0]
                for y in ys:
                    self.fig.add_hline(y=y, line_color="red", row=i+1, col=1)

        print("A", time.time()-t0)
        t0 = time.time()
        

        row = len(self.tracks)+1
        for i,stat in enumerate(layer_stats):
            self.fig.update_yaxes(title_text=stat, row=row, col=1)
            for j,track in enumerate(self.tracks.aln_tracks):
                stats = self.tracks.refstats[track.name,self.prms.layer,stat]
                self.fig.add_trace(go.Scattergl(
                    name=track.name,
                    legendgroup=track.name,
                    showlegend=i==0,
                    x=stats.index,
                    y=stats,
                    line={"color":self.prms.track_colors[j]},
                ), row=row, col=1)
            row += 1

        print("B", time.time()-t0)
        t0 = time.time()

        for stat in cmp_stats:
            self.fig.update_yaxes(title_text=stat + " stat", row=row, col=1)
            stats = self.tracks.refstats[stat,self.prms.layer,"stat"]
            self.fig.add_trace(go.Scattergl(
                name=stat,
                x=stats.index,
                y=stats,
                line={"color":"red"},
            ), row=row, col=1)
            row += 1

        print("C", time.time()-t0)
        t0 = time.time()

        if self.prms.select_ref is not None:
            self.fig.add_vline(x=self.prms.select_ref, line_color="red")

        cax = {"colorbar" : {
            "title" : layer_label,
             "len" : 250, "y" : 1,
             "lenmode" : "pixels",
             "yanchor" : "top"}}
        cax.update(LAYER_COLORS[self.prms.layer])

        self.fig.update_xaxes(side='top', showticklabels=True, row=1, col=1)
        self.fig.update_xaxes(showticklabels=True, row=n_rows, col=1)

        height = max(500, 100*np.sum(row_heights))

        self.fig.update_layout(
            coloraxis=cax, dragmode="pan", 
            autosize=True, height=height,
            margin={"t":50},
            legend={"x":1,"y":0,"xanchor":"left"}
        )
        self.fig.update_layout()

        self.fig.update_yaxes(showticklabels=False)#, title="Reads")
        print("D", time.time()-t0)

    def show(self):
        fig_conf = {"scrollZoom" : True, "displayModeBar" : True}

        if self.prms.outfile is not None:
            self.fig.write_html(self.prms.outfile, config=fig_conf)
        else:
            self.fig.show(config=fig_conf)

OPTS = (
    Opt("ref_bounds", "track_io", type=str_to_coord),
    Opt("input", "track_io", nargs="+"),
    Opt("layer", "trackplot"),
    Opt(("-r", "--refstats"), "track_io", default="mean", type=comma_split),
    Opt(("-f", "--full-overlap"), "track_io", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
)

def main(conf):
    conf.track_io.layers.append(conf.trackplot.layer)
    conf.track_io.refstats_layers = [conf.trackplot.layer]
    Trackplot(conf=conf).show()
