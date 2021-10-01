import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from ...dtw.track import LAYER_META
from ...index import str_to_coord
from ...dtw.track_io import TrackIO
from ...argparse import Opt, comma_split

def trackplot(tracks, layer):

    names = [t.name for t in tracks]
    fig = make_subplots(
        rows=len(tracks), cols=1, 
        subplot_titles=names, 
        shared_xaxes=True, 
        x_title=tracks[0].coords.ref_name,#)
        #y_title="Reads",
        vertical_spacing=0.125/len(tracks))

    layer_desc = LAYER_META[layer].label

    for i,track in enumerate(tracks):
        mat = track.mat["dtw"][layer]
        hover = "<br>".join([
            track.coords.ref_name + ":%{x}",
            "Read: %{y}",
            layer_desc + ": %{z}"])

        fig.add_trace(go.Heatmap(
            name=track.name,
            x=mat.columns,
            y=track.alignments["read_id"],
            z=mat,
            zsmooth=False,
            hovertemplate=hover,
            coloraxis="coloraxis",
        ), row=i+1, col=1)

    fig.update_layout(
        coloraxis={
            "colorscale": "balance", 
            "cmid" : 0,
            "colorbar": {"title" : layer_desc,
                         "len" : 250,
                         "lenmode" : "pixels",
                         "y" : 1,
                         "yanchor" : "top"}},
        autosize=False, height=max(500, 300*len(tracks)), width=800
    )

    fig.update_yaxes(showticklabels=False)#, title="Reads")
    return fig

OPTS = (
    Opt("ref_bounds", "track_io", type=str_to_coord),
    Opt("input", "track_io", nargs="+"),
    Opt(("-f", "--full-overlap"), "track_io", action="store_true"),
    Opt("layer", "track_io", "layers", type=str, default="current"),
)

def main(conf):
    io = TrackIO(conf=conf)

    tracks = io.load_refs(load_mat=True)
    fig = trackplot(tracks, conf.track_io.layers)
    fig.write_html('first_figure.html')#, auto_open=True)
