import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from ... import nt

from ...dtw.track import LAYER_META
from ...index import str_to_coord
from ...dtw.track_io import TrackIO
from ...argparse import Opt, comma_split
from ...fast5 import parse_read_ids
from ...sigproc import ProcRead

def dotplot(conf, read, tracks):
    fig = make_subplots(
        rows=2, cols=1, 
        row_heights=[1,3],
        vertical_spacing=0.01,
        shared_xaxes=True)

    for i,track in enumerate(tracks):
        for aln_id, aln in track.alignments.iterrows():
            dtw = track.layers["dtw"].loc[aln_id]

            samp_min = dtw["start"].min()
            max_i = dtw["start"].argmax()
            samp_max = dtw["start"].iloc[max_i] + dtw["length"].iloc[max_i]
            samps = np.arange(samp_min, samp_max)
            raw_norm = read.get_norm_signal(samp_min, samp_max)
            mask = ((raw_norm >= 40) &
                    (raw_norm <= conf.event_detector.max_mean))

            kmers = track.coords.kmers[dtw.index]
            model_current = track.model[kmers]
            print(track.conf.pore_model.name)

            aln_bases = nt.kmer_base(kmers, 2)

            ymin = min(np.min(model_current), np.min(raw_norm[mask]))
            ymax = max(np.max(model_current), np.max(raw_norm[mask]))

            base_colors = ["#80ff80", "#8080ff", "#ffbd00", "#ff8080"] #A,C,G,T
            for base, color in enumerate(base_colors):
                base_dtw = dtw[aln_bases == base]
                starts = base_dtw['start']
                ends = starts + base_dtw['length'] - 1
                nones = [None]*len(base_dtw)

                ys = [ymax,ymax,ymin,ymin,None]*len(base_dtw)
                xs = np.dstack([starts, ends, ends, starts, nones]).reshape(len(ys))
                
                fig.add_trace(go.Scatter(
                    x=xs,y=ys, fill="toself",
                    fillcolor=color,
                    hoverinfo="skip",
                    line={"width" : 0},
                    name=nt.base_to_char(base)
                ), row=1, col=1)

            fig.add_trace(go.Scattergl(
                x=samps[mask], y=raw_norm[mask],
                hoverinfo="skip",
                name="Raw Signal",
                mode="markers",
                marker={"size":2, "color":"black"}
            ), row=1, col=1)

            fig.add_trace(go.Scattergl(
                name = "Model Current",
                x=dtw["start"], y=model_current,
                line={"color":"white", "width":2, "shape" : "hv"}
            ), row=1, col=1)

            fig.add_trace(go.Scatter(
                x=dtw["start"], y=dtw.index,
                name=track.name,
                line={"color":"purple", "width":2, "shape" : "hv"}
            ), row=2, col=1)

    fig.update_layout(yaxis={"fixedrange" : True}, dragmode="pan")#, scroll_zoom=True)

    #fig.update_layout(
    #    coloraxis={
    #        "colorscale": "balance", 
    #        "cmid" : 0,
    #        "colorbar": {"title" : layer_desc,
    #                     "len" : 250,
    #                     "lenmode" : "pixels",
    #                     "y" : 1,
    #                     "yanchor" : "top"}},
    #    autosize=False, height=max(500, 300*len(tracks)), width=800
    #)

    #fig.update_yaxes(showticklabels=False)#, title="Reads")

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
    fast5s = io.get_fast5_reader()

    for read_id, tracks in io.iter_reads():
        fast5_read = fast5s[read_id]
        if isinstance(fast5_read, ProcRead):
            read = fast5_read
        else:
            read = ProcRead(fast5_read, conf=io.conf)

        fig = dotplot(conf, read, tracks)

        #fig.show()

        fig.write_html(conf.out_prefix + read_id + ".html", config={"scrollZoom" : True, "displayModeBar" : True})
