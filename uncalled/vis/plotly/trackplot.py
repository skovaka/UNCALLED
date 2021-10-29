import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time

from ... import config
from ...dtw.aln_track import LAYERS, parse_layer
from ...index import str_to_coord
from ...dtw.tracks import Tracks, REFSTAT_LABELS
from ...argparse import Opt, comma_split

class TrackplotParams(config.ParamGroup):
    _name = "trackplot"
TrackplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("layer", "current", None, "Layer to plot"),
    ("track_colors", ["#AA0DFE", "#1CA71C", "#4676FF", "red"], list, ""),
    ("select_ref", None, str, "Reference Selection"),
    ("select_read", None, str, "Read Selection"),
    ("width", None, int, "Figure width"),
    ("min_height", 700, int, "Minimum figure height"),
    ("track_height", 100, int, "Minimum per-track figure height"),
    ("outfile", None, str, "Output file"),
    #("track_colors", ["#AA0DFE", "#1CA71C", "#6A76FC"], list, ""),
)

LAYER_COLORS = {
    ("dtw", "model_diff") : {"colorscale" : "RdBu", "cmid" : 0, "cmax" : 20, "cmin" : -20, "reversescale":True},
    ("dtw", "current") : {"colorscale" : "viridis"},
    ("dtw", "dwell") : {"colorscale" : "viridis", "cmin" : 0, "cmax" : 25},
    ("cmp", "mean_ref_dist") : {"colorscale" : "RdYlGn", "cmin" : 0, "cmid" :5, "cmax" : 10, "reversescale":True},
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
            cmp_stats = []
        else:
            refstat_tracks = self.tracks.refstats.columns.get_level_values("track").unique()
            cmp_stats = refstat_tracks.difference(names)

            #layer_stats = ["box"]
            if len(refstat_tracks.intersection(names)) == len(names):
                layer_stats = self.tracks.refstats[names].columns.get_level_values("stat").unique()
                print(self.tracks.refstats)
            else:
                layer_stats = pd.Index([])

        n_stats = len(layer_stats)+len(cmp_stats)
        row_heights = [1]*n_stats + [2]*len(self.tracks)
        n_rows = len(row_heights)

        t0 = time.time()
        ref_title = "Reference (%s)" % self.tracks.coords.ref_name
        self.fig = make_subplots(
            rows=n_rows, cols=1, 
            #subplot_titles=[ref_title], 
            row_heights=row_heights,
            shared_xaxes=True, 
            x_title=ref_title,
            #y_title="Reads",
            vertical_spacing=0.125/n_rows)

        group,layer = self.prms.layer
        layer_label = LAYERS[group][layer].label

        t0 = time.time()
        for i,track in enumerate(self.tracks.all):
            mat = track.mat[self.prms.layer]#.dropna(how="all")
            hover = "<br>".join([
                track.coords.ref_name + ":%{x}",
                #"Read: %{y}",
                layer_label + ": %{z}"])

            row = n_stats+i+1
            self.fig.add_trace(go.Heatmap(
                name=track.desc,
                x=mat.columns,
                #y=track.alignments["read_id"],
                z=mat,
                zsmooth=False,
                hovertemplate=hover,
                coloraxis="coloraxis",
            ), row=row, col=1)

            self.fig.update_yaxes(
                title_text=track.desc, 
                #title_font_color=self.prms.track_colors[i], 
                showticklabels=False,
                row=row, col=1)

            if self.prms.select_read is not None:
                ys = np.where((self.prms.select_read == track.alignments["read_id"]).to_numpy())[0]
                for y in ys:
                    self.fig.add_hline(y=y, line_color="red", row=row, col=1)

        t0 = time.time()
        

        row = 1
        for i,stat in enumerate(layer_stats):

            if stat == "box":
                self.fig.update_yaxes(title_text=layer_label, row=row, col=1)
                for j,track in enumerate(self.tracks.aln_tracks):
                    stats = self.tracks.refstats[track.name,group,layer]
                    print(stats)
                    for idx in stats.index[:-1]:
                        self.fig.add_vline(x=idx-0.5, line_color="black", row=row, col=1)
                    self.fig.add_trace(go.Box(
                        x=stats.index - 0.25 + j*0.5,
                        median=stats["median"],
                        lowerfence=stats["q5"],
                        q1=stats["q25"],
                        q3=stats["q75"],
                        upperfence=stats["q95"],
                        #sd=stats["stdv"],
                        name=track.desc,
                        #fillcolor=self.prms.track_colors[j],
                        line_color=self.prms.track_colors[j]
                    ), row=row, col=1)
                row += 1
                continue

            print(i,stat)
            self.fig.update_yaxes(title_text=REFSTAT_LABELS[stat], row=row, col=1)
            for j,track in enumerate(self.tracks.aln_tracks):
                group,layer = self.prms.layer
                stats = self.tracks.refstats[track.name,group,layer,stat]
                self.fig.add_trace(go.Scattergl(
                    name=track.desc,
                    legendgroup=track.desc,
                    showlegend=i==0,
                    x=stats.index,
                    y=stats,
                    opacity=0.5,
                    mode="markers", marker={"size":3,"color":self.prms.track_colors[j]},
                    #line={"color":self.prms.track_colors[j]},
                ), row=row, col=1)
            row += 1

        t0 = time.time()

        for stat in cmp_stats:
            self.fig.update_yaxes(title_text=REFSTAT_LABELS[stat] + " Stat", row=row, col=1)
            group,layer = self.prms.layer
            stats = self.tracks.refstats[stat,group,layer,"stat"]
            self.fig.add_trace(go.Scattergl(
                name="Compare",
                x=stats.index,
                y=stats,
                mode="lines",
                line={"color":"red"},
            ), row=row, col=1)
            row += 1

        t0 = time.time()

        if self.prms.select_ref is not None:
            self.fig.add_vline(x=self.prms.select_ref, line_color="red")

        cax = {"colorbar" : {
            "title" : layer_label,
             "y" : 0.5, "len" : 0.5, 
             "yanchor" : "top"}}

        if self.prms.layer in LAYER_COLORS:
            cax.update(LAYER_COLORS[self.prms.layer])

        self.fig.update_xaxes(side='top', showticklabels=True, row=1, col=1)
        self.fig.update_xaxes(showticklabels=True, row=n_rows, col=1)

        height = max(700, 100*np.sum(row_heights))

        self.fig.update_layout(
            coloraxis=cax, dragmode="pan", 
            autosize=True, height=height,
            margin={"t":50},
            legend={"x":1,"y":1,"xanchor":"left"}
        )
        self.fig.update_layout()

    def show(self):
        fig_conf = {
            "toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None, "scale" : 2},
            "scrollZoom" : True, 
            "displayModeBar" : True}

        if self.prms.outfile is not None:
            self.fig.write_html(self.prms.outfile, config=fig_conf)
        else:
            self.fig.show(config=fig_conf)

OPTS = (
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt("input", "tracks", nargs="+"),
    Opt("layer", "trackplot", type=parse_layer),
    Opt(("-r", "--refstats"), "tracks", default=None, type=comma_split),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
)

def main(conf):
    conf.tracks.layers.append(conf.trackplot.layer)
    conf.tracks.refstats_layers = [conf.trackplot.layer]
    Trackplot(conf=conf).show()
