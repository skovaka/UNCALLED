import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time

from ... import config
from ...dtw.track import LAYERS
from ...index import str_to_coord
from ...dtw.track_io import TrackIO
from ...argparse import Opt, comma_split

class TrackplotParams(config.ParamGroup):
    _name = "trackplot"
TrackplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("layer", "current", str, "Layer to plot"),
    ("select_ref", None, str, "Reference Selection"),
    ("select_read", None, str, "Read Selection"),
    ("outfile", None, str, "Output file"),
    #("track_colors", ["#AA0DFE", "#1CA71C", "#6A76FC"], list, ""),
)

LAYER_COLORS = {
    "model_diff" : {"colorscale" : "balance", "cmid" : 0},
    "current" : {"colorscale" : "viridis"},
    "dwell" : {"colorscale" : "viridis", "cmin" : 0, "cmax" : 50},
}

class Trackplot:

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("trackplot", *args, **kwargs)

        if self.prms.tracks is None:
            self.tracks = TrackIO(conf=self.conf)
            self.tracks.load_refs(load_mat=True)
        else:
            self.tracks = self.prms.tracks

        t0 = time.time()
        names = [t.name for t in self.tracks]
        self.fig = make_subplots(
            rows=len(self.tracks), cols=1, 
            subplot_titles=names, 
            shared_xaxes=True, 
            x_title=self.tracks.coords.ref_name,#)
            #y_title="Reads",
            vertical_spacing=0.125/len(self.tracks))

        layer_label = LAYERS["dtw"][self.prms.layer].label

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

            if self.prms.select_read is not None:
                ys = np.where((self.prms.select_read == track.alignments["read_id"]).to_numpy())[0]
                for y in ys:
                    self.fig.add_hline(y=y, line_color="red", row=i+1, col=1)

        if self.prms.select_ref is not None:
            self.fig.add_vline(x=self.prms.select_ref, line_color="red")

        cax = {"colorbar" : {
            "title" : layer_label,
             "len" : 250, "y" : 1,
             "lenmode" : "pixels",
             "yanchor" : "top"}}
        cax.update(LAYER_COLORS[self.prms.layer])

        self.fig.update_layout(
            coloraxis=cax, dragmode="pan", 
            autosize=True, height=max(500, 300*len(self.tracks)),
            margin={"t":50},
            #paper_bgcolor="lightgray",
        )

        self.fig.update_yaxes(showticklabels=False)#, title="Reads")

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
    Opt("select_read", "trackplot"),
    Opt(("-f", "--full-overlap"), "track_io", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
)

def main(conf):
    conf.track_io.layers.append(conf.trackplot.layer)
    Trackplot(conf=conf).show()
