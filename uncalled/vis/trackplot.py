import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time
import sys

from .. import config
from ..dtw.aln_track import LAYERS, parse_layer, parse_layers
from ..index import str_to_coord
from ..dtw.tracks import Tracks, REFSTAT_LABELS, COMPARE_REFSTATS
from ..argparse import Opt, comma_split

class TrackplotParams(config.ParamGroup):
    _name = "trackplot"
TrackplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("panels", None, None, "List of tuples specifying which panels to display. First element of each tuple specifies plot type (e.g. mat, box, line, etc.), second element specifies which layer(s) to display."),
    ("track_colors", ["#AA0DFE", "#1CA71C", "#4676FF", "red"], list, ""),
    ("select_ref", None, str, "Reference Selection"),
    ("select_read", None, str, "Read Selection"),
    ("width", None, int, "Figure width"),
    ("panel_heights", None, None, "Relative height of each panel"),
    ("min_height", 700, int, "Minimum figure height"),
    ("track_height", 100, int, "Minimum per-track figure height"),
    ("outfile", None, str, "Output file"),
)

_CMP_COLOR = {"colorscale" : "RdYlGn", "cmin" : 0, "cmid" :5, "cmax" : 10, "reversescale":True}

LAYER_COLORS = {
    ("dtw", "model_diff") : {"colorscale" : "RdBu", "cmid" : 0, "cmax" : 20, "cmin" : -20, "reversescale":True},
    ("dtw", "current") : {"colorscale" : "viridis"},
    ("dtw", "dwell") : {"colorscale" : "viridis", "cmin" : 0, "cmax" : 25},
    ("cmp", "mean_ref_dist") : _CMP_COLOR,
    ("bc_cmp", "mean_ref_dist") : _CMP_COLOR,
}

LAYER_PANELS = {"mat", "box"}
REFSTAT_PANELS = {"line", "scatter"}

MULTIROW_PANEL = {
    "mat" : True, "box" : False, "line" : False, "scatter" : False
}

DEFAULT_HEIGHTS = {
    "mat" : 3, "box" : 2, "line" : 1, "scatter" : 1
}

class Trackplot:

    def _parse_panels(self):
        layers = list()
        ref_layers = list()
        ref_stats = set()
        for panel, layer in self.prms.panels:
            if panel in LAYER_PANELS:
                layers.append(layer)
                if panel == "box":
                    ref_layers.append(layer)
                    ref_stats.update(["median","stdv","q5","q25","q75","q95"])
            elif panel in REFSTAT_PANELS:
                spl = layer.split(".")
                ref_layers.append(".".join(spl[:-1]))
                ref_stats.add(spl[-1])
        
        prms = self.conf.tracks
        prms.layers = list(parse_layers(layers))
        prms.refstats_layers = list(parse_layers(ref_layers))
        prms.refstats = list(ref_stats)

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("trackplot", *args, **kwargs)

        self._parse_panels()

        if self.prms.tracks is None:
            self.tracks = Tracks(conf=self.conf)
            self.tracks.load_refs(load_mat=True)
        else:
            self.tracks = self.prms.tracks

        names = [t.name for t in self.tracks]

        if self.prms.panel_heights is None:
            panel_heights = list()
            for panel,_ in self.prms.panels:
                if MULTIROW_PANEL[panel]:
                    panel_heights += [DEFAULT_HEIGHTS[panel]] * len(self.tracks)
                else:
                    panel_heights.append(DEFAULT_HEIGHTS[panel])
        else:
            panel_heights = self.prms.panel_heights

        n_rows = len(panel_heights)

        t0 = time.time()
        ref_title = "Reference (%s)" % self.tracks.coords.ref_name
        self.fig = make_subplots(
            rows=n_rows, cols=1, 
            row_heights=self.prms.panel_heights,
            shared_xaxes=True, 
            x_title=ref_title,
            #y_title="Reads",
            vertical_spacing=0.125/n_rows)

        row = 1
        for panel,layer in self.prms.panels:
            fn = getattr(self, "_"+panel)
            fn(row, layer)
            row += len(self.tracks) if MULTIROW_PANEL[panel] else 1

        if self.prms.select_ref is not None:
            self.fig.add_vline(x=self.prms.select_ref, line_color="red")

        #cax = {"colorbar" : {
        #    "title" : layer_label,
        #     "y" : 0.5, "len" : 0.5, 
        #     "yanchor" : "top"}}
        #if self.prms.layer in LAYER_COLORS:
        #    cax.update(LAYER_COLORS[self.prms.layer])
        #self.fig.update_layout(coloraxis=cax)

        self.fig.update_xaxes(side='top', showticklabels=True, row=1, col=1)
        self.fig.update_xaxes(showticklabels=True, row=n_rows, col=1)

        height = max(700, 100*np.sum(panel_heights))

        self.fig.update_layout(
            dragmode="pan", 
            autosize=True, height=height,
            margin={"t":50},
            legend={"x":1,"y":1,"xanchor":"left"}
        )
        
    def _mat(self, row, layer):
        group,layer = parse_layer(layer)

        layer_label = LAYERS[group][layer].label

        t0 = time.time()
        for i,track in enumerate(self.tracks.aln_tracks):
            self.fig.update_yaxes(
                title_text=f"{track.desc} ({layer_label})", 
                showticklabels=False,
                row=row, col=1)

            if track.empty: 
                sys.stderr.write("Warning: no alignments loaded from track \"%s\"\n" % track.desc)
                row += 1
                continue

            mat = track.mat[(group,layer)]#.dropna(how="all")
            hover = "<br>".join([
                track.coords.ref_name + ":%{x}",
                #"Read: %{y}",
                layer_label + ": %{z}"])

            self.fig.add_trace(go.Heatmap(
                name=track.desc,
                x=mat.columns,
                #y=track.alignments["read_id"],
                z=mat,
                zsmooth=False,
                hovertemplate=hover,
                coloraxis="coloraxis",
            ), row=row, col=1)

            if self.prms.select_read is not None:
                ys = np.where((self.prms.select_read == track.alignments["read_id"]).to_numpy())[0]
                for y in ys:
                    self.fig.add_hline(y=y, line_color="red", row=row, col=1)
            row += 1

        cax = {"colorbar" : {
            "title" : layer_label,
             "y" : 0.5, "len" : 0.5, 
             "yanchor" : "top"}}
        if (group,layer) in LAYER_COLORS:
            cax.update(LAYER_COLORS[(group,layer)])
        self.fig.update_layout(coloraxis=cax)
        
    def _box(self, row, layer):
        group,layer = parse_layer(layer)

        layer_label = LAYERS[group][layer].label

        self.fig.update_yaxes(title_text=layer_label, row=row, col=1)
        for j,track in enumerate(self.tracks.aln_tracks):
            stats = self.tracks.refstats[track.name,group,layer]
            for idx in stats.index[:-1]:
                self.fig.add_vline(x=idx-0.5, line_color="black", row=row, col=1)
            self.fig.add_trace(go.Box(
                x=stats.index - 0.25 + j*0.5,
                median=stats["median"],
                lowerfence=stats["q5"],
                q1=stats["q25"],
                q3=stats["q75"],
                upperfence=stats["q95"],
                name=track.desc,
                #fillcolor=self.prms.track_colors[j],
                line_color=self.prms.track_colors[j]
            ), row=row, col=1)

    def _line(self, row, layer):
        self._refstat(row, layer, "line")#line_color=self.prms.track_colors[j])

    def _scatter(self, row, layer):
        self._refstat(row, layer, "scatter")#mode="markers", marker={"size" : 3, "color" : self.prms.track_colors[j]})

    def _stat_kw(self, plot, color):
        if plot == "line":
            return {"line_color" : color}
        elif plot == "scatter":
            return {"mode" : "markers", 
                  "marker" : {
                    "size" : 3, "opacity" : 0.5, 
                    "color" : color}}
        raise ValueError( f"Unknown plot type: {plot}")

    def _refstat(self, row, layer, plot):
        spl = layer.split(".")
        group, layer = parse_layer(".".join(spl[:-1]))
        stat = spl[-1]

        layer_label = LAYERS[group][layer].label
        stat_label = REFSTAT_LABELS[stat]

        if stat in COMPARE_REFSTATS:
            self.fig.update_yaxes(title_text=f"{layer_label} {stat_label}", row=row, col=1)
            stats = self.tracks.refstats[stat,group,layer,"stat"]
            self.fig.add_trace(go.Scattergl(
                name="Compare",
                x=stats.index,
                y=stats,
                **self._stat_kw(plot, "red")
            ), row=row, col=1)
            return

        self.fig.update_yaxes(title_text=REFSTAT_LABELS[stat], row=row, col=1)
        for j,track in enumerate(self.tracks.aln_tracks):
            stats = self.tracks.refstats[track.name,group,layer,stat]

            self.fig.add_trace(go.Scattergl(
                name=track.desc,
                legendgroup=track.desc,
                #showlegend=i==0,
                x=stats.index,
                y=stats,
                **self._stat_kw(plot, self.prms.track_colors[j])
            ), row=row, col=1)

        #for stat in cmp_stats:

    def show(self):
        fig_conf = {
            "toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None, "scale" : 2},
            "scrollZoom" : True, 
            "displayModeBar" : True}

        if self.prms.outfile is not None:
            self.fig.write_html(self.prms.outfile, config=fig_conf)
        else:
            self.fig.show(config=fig_conf)

def mat_opt(layer):
    return ("mat", layer)
    
def panel_opt(name):
    return (lambda arg: (name, arg))

OPTS = (
    Opt("input", "tracks", nargs="+"),
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-H", "--panel-heights"), "trackplot", nargs="+", type=int),

    Opt("--mat", dest="panels",
        metavar="LAYER", action="append", type=panel_opt("mat"),
        help="Display a ref-by-read matrix of specified alignment layer"), 

    Opt("--box", dest="panels", #"trackplot", "panels", 
        metavar="LAYER", action="append", type=panel_opt("box"),
        help="Display a boxplot of specified layer"), 

    Opt("--line", dest="panels", #"trackplot", "panels", 
        metavar="LAYER.STAT", action="append", type=panel_opt("line"),
        help="Display a line plot of specifed layer summary statistic"), 

    Opt("--scatter", dest="panels", #"trackplot", "panels", 
        metavar="LAYER.STAT", action="append", type=panel_opt("scatter"),
        help="Display a line plot of specifed layer summary statistic"), 
    Opt(("-o", "--outfile"), "trackplot"),
)

def main(conf):
    """Plot alignment tracks and per-reference statistics"""
    #conf.tracks.layers.append(conf.trackplot.layer)
    #conf.tracks.refstats_layers = [conf.trackplot.layer]
    conf.trackplot.panels = conf.panels
    Trackplot(conf=conf).show()
