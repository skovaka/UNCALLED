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

#cmap = px.colors.qualitative.Dark24[1:]
cmap = ["#AA0DFE", "#1CA71C"]#, "#6A76FC"]
track_colors = cmap

def dotplot(conf, tracks):
    fig = make_subplots(
        rows=2, cols=2, 
        row_heights=[1,3],
        column_widths=[3,1],
        vertical_spacing=0.01,
        horizontal_spacing=0.01,
        shared_xaxes=True,
        shared_yaxes=True)

    Sigplot(tracks, track_colors=track_colors, conf=conf).plot(fig)

    hover_df = dict()

    for i,track in enumerate(tracks):

        track_df = list()
        
        #TODO write AlnIO AlnTracks.tracks
        for aln_id, aln in track.alignments.iterrows():
            dtw = track.layers["dtw"].xs(aln_id, level="aln_id")

            fwd, refs = track.coords.mref_to_ref(dtw.index)

            track_df.append(pd.DataFrame({
                "ref" : refs,
                "sample_mid" : dtw["start"] + dtw["length"]/2,
                "current" : dtw["current"],
                "dwell" : dtw["dwell"],
                "model_diff" : dtw["model_diff"],
            }).set_index("ref", drop=True))

            fig.add_trace(go.Scatter(
                x=dtw["start"], y=refs,
                name=track.name,
                legendgroup=track.desc,
                line={"color":track_colors[i], "width":2, "shape" : "hv"},
                hoverinfo="skip",
            ), row=2, col=1)

            fig.add_trace(go.Scatter(
                x=dtw["model_diff"], y=refs-0.5, 
                name=track.name, legendgroup=track.desc,
                line={"color":track_colors[i], "width":2, "shape" : "hv"},
            ), row=2, col=2)

        hover_df[track.name] = pd.concat(track_df)

    hover_df = pd.concat(hover_df, axis=1)
    hover_coords = hover_df.xs("sample_mid", axis=1, level=1).mean(axis=1)
    hoverdata = hover_df.drop("sample_mid", axis=1, level=1).to_numpy()

    hover_rows = [
        track.coords.ref_name + ":%{y:,d}"
    ]
    labels = [
        "Current (pA): ", 
        "Dwell (ms): ", 
        "Model pA Diff: "]
    for i,label in enumerate(labels):
        s = "<b>"+label+"</b>"
        #s = "<b style=\"font-family:mono\">"+label+"</b>"
        fields = list()
        for j in range(len(tracks)):
            fields.append(
                '<span style="color:%s;float:right">%%{customdata[%d]:.2f}</span>' % 
                (track_colors[j], j*3+i))
        hover_rows.append(s + ", ".join(fields))

    #hovertemplate="<br>".join([
    #    "<b>" + track.coords.ref_name + ":</b> %{customdata[0]:d}",
    #    "<b>Raw Sample:</b> <span style=\"color:red\">%{customdata[1]:d}</span>",
    #    "<b>Current (pA):</b> %{customdata[2]:.2f}",
    #    "<b>Dwell Time (pA):</b> %{customdata[3]:.2f}",
    #    "<b>Model pA Difference:</b> %{customdata[4]:.2f}"
    #])

    fig.add_trace(go.Scatter(
        x=hover_coords, y=hover_coords.index,
        mode="markers", marker={"size":0,"color":"rgba(0,0,0,0)"},
        name="",
        customdata=hoverdata,
        hovertemplate="<b>"+"<br>".join(hover_rows)+"</b>",
        hoverlabel={"bgcolor":"rgba(255,255,255,1)"},
        showlegend=False
    ), row=2,col=1)
        #customdata=hoverdata.to_numpy(),


    if fwd == conf.is_rna:
        fig.update_yaxes(autorange="reversed", row=2, col=1)
        fig.update_yaxes(autorange="reversed", row=2, col=2)

    fig.update_yaxes(
        title_text="Reference (%s)" % aln["ref_name"], 
        #showspikes=True,
        row=2, col=1)

    axis_kw = dict(
        showspikes=True,
        spikemode="across",
        spikecolor="darkgray",
        spikethickness=1)

    fig.update_xaxes(**axis_kw)
    fig.update_yaxes(**axis_kw)
    #fig.update_yaxes(showspikes=True)

    fig.update_xaxes(
        title_text="Raw Sample", 
        tickformat="d", 
        #showspikes=True,
        row=2, col=1)

    fig.update_layout(
        #hovermode="x unified",
        hoverdistance=20,
        dragmode="pan", 
        legend={"bgcolor" : "#e6edf6"})#, scroll_zoom=True)

    #fig.update_traces(xaxis='x2')

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
