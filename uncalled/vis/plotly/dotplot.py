import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from .sigplot import Sigplot

from ... import nt

from ...dtw.track import LAYER_META
from ...index import str_to_coord
from ...dtw.track_io import TrackIO
from ...argparse import Opt, comma_split
from ...fast5 import parse_read_ids
from ...sigproc import ProcRead

def dotplot(conf, tracks):
    fig = make_subplots(
        rows=2, cols=1, 
        row_heights=[1,3],
        vertical_spacing=0.01,
        shared_xaxes=True)

    Sigplot(tracks, conf=conf).plot(fig)

    for i,track in enumerate(tracks):
        for aln_id, aln in track.alignments.iterrows():

            dtw = track.layers["dtw"].loc[aln_id]

            fig.add_trace(go.Scatter(
                x=dtw["start"], y=dtw.index,
                name=track.name,
                line={"color":"purple", "width":2, "shape" : "hv"}
            ), row=2, col=1)

    fig.update_layout(yaxis={"fixedrange" : True}, dragmode="pan")#, scroll_zoom=True)

    return fig

OPTS = (
    Opt("input", "track_io", nargs="+"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "track_io", type=str_to_coord),
    Opt(("-l", "--read-filter"), "track_io", type=parse_read_ids),
    #Opt(("-L", "--layers"), "dotplot", type=comma_split),
)

def main(conf):
    """plot a dotplot"""
    io = TrackIO(conf=conf)

    for read_id, tracks in io.iter_reads():
        print(read_id)
        fig = dotplot(io.conf, tracks)

        fig.write_html(conf.out_prefix + read_id + ".html", config={"scrollZoom" : True, "displayModeBar" : True})
