import plotly.express as px
import pandas as pd
import numpy as np

import dash
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State

from .trackplot import Trackplot
from ... import config
from ...index import str_to_coord
from ...dtw.track_io import TrackIO
from ...dtw.track import LAYERS
from ...argparse import Opt, comma_split


OPTS = (
    Opt("ref_bounds", "track_io", type=str_to_coord),
    Opt("input", "track_io", nargs="+"),
    Opt("layer", "trackplot", default="model_diff", nargs="?"),
    Opt(("-f", "--full-overlap"), "track_io", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
)

def main(conf):
    """Nanopore DTW genome browser"""
    tracks = TrackIO(conf=conf)
    print("Loading tracks...")
    tracks.load_refs(load_mat=True)
    print("Done")
    browser(tracks, conf)

def browser(tracks, conf):
    external_stylesheets = ["https://www.w3schools.com/w3css/4/w3.css"]

    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
    app.title = "Uncalled4 Browser"

    app.layout = html.Div(children=[
        html.Div(
            html.H2(html.B("Uncalled4 Genome Browser")), 
            className="w3-container w3-deep-purple"),

        html.Div([
            html.Div([
                dcc.Graph(
                    id="trackplot",
                    config = {"scrollZoom" : True, "displayModeBar" : True}
                )],
            className="w3-container w3-twothird w3-card"),

            html.Div([
                html.B("Active layer: "), 
                dcc.Dropdown(
                    options=[
                        {"label": "Current (pA)", "value": "current"},
                        {"label": "Dwell Time (ms)", "value": "dwell"},
                        {"label": "Model pA Difference", "value": "model_diff"},
                    ],
                    value=conf.trackplot.layer, 
                    clearable=False, multi=False,
                    id="trackplot-layer"),
                html.Div([], id="click-info", className="w3-container w3-card"),
            ], className="w3-container w3-third")
        ])
    ])

    #@app.callback(
    #    Output("trackplot", "figure"),
    #    Input("trackplot-layer", "value"))
    #def update_trackplot(layer):
    #    fig = Trackplot(tracks, layer, conf=conf).fig
    #    return fig

    @app.callback(
        Output("trackplot", "figure"),
        Output("click-info", "children"),
        Input("trackplot-layer", "value"),
        Input("trackplot", "clickData"))
    def update_trackplot(layer, click):
        if click is None:
            ref = aln = read = None
            info = ""
        else:
            coord = click["points"][0]
            ref = coord["x"]
            track = tracks.all[coord["curveNumber"]]
            aln = track.alignments.iloc[coord["y"]]
            mref = track.coords.ref_to_mref(ref, aln.fwd)
            read = aln["read_id"]

            layers = track.layers.loc[(mref, aln.name)]["dtw"]

            table = list()
            table.append(html.Tr(html.Td(html.B("%s:%d" % (tracks.coords.ref_name, ref)), colSpan=2)))
            table.append(html.Tr(html.Td([html.B("Read "), read], colSpan=2)))
            for l in ["current", "dwell", "model_diff"]:
                table.append(html.Tr([
                    html.Td(html.B(LAYERS["dtw"][l].label)), 
                    html.Td("%.3f"%layers[l], style={"text-align":"right"})]))

            info = [html.Table(table)]
            
        fig = Trackplot(
            tracks, layer, 
            select_ref=ref, select_read=read, 
            conf=conf).fig
        fig.update_layout(uirevision=True)

        return fig, info

    app.run_server(debug=True)
