import plotly.express as px
import pandas as pd
import numpy as np

import dash
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State

from .trackplot import Trackplot
from .. import config
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..dtw.aln_track import LAYERS, parse_layer
from ..argparse import Opt, comma_split


OPTS = (
    Opt("input", "tracks", nargs="+"),
    Opt("ref_bounds", "tracks", type=str_to_coord),
    #Opt("layer", "trackplot", default="current", nargs="?"),
    Opt(("-r", "--refstats"), "tracks", default=None, type=comma_split),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
)

def main(conf):
    """Nanopore DTW genome browser"""
    conf.tracks.layers.append("cmp.mean_ref_dist")
    conf.tracks.refstats_layers.append("cmp.mean_ref_dist")
    tracks = Tracks(conf=conf)
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
                dcc.Graph(#[dcc.Loading(type="circle"),
                    id="trackplot",
                    config = {"scrollZoom" : True, "displayModeBar" : True}
            )], className="w3-container w3-twothird w3-card"),

            html.Div([
                html.Div([
                    html.B("Active layer: "), 
                    dcc.Dropdown(
                        options=[
                            {"label": "Current (pA)", "value": "current"},
                            {"label": "Dwell Time (ms)", "value": "dwell"},
                            {"label": "Model pA Difference", "value": "model_diff"},
                            {"label": "Mean Ref. Dist", "value": "cmp.mean_ref_dist"},
                        ],
                        value=conf.trackplot.layer, 
                        clearable=False, multi=False,
                        id="trackplot-layer"),
                ]),

                html.Div([
                    html.Table([], id="info-table"),
                    html.Button("View Dotplot", id="dotplot-btn", style={"margin" : "5px"})
                    ], style={"display" : "none"},
                    id="selection-card", className="w3-container w3-card"),

            ], className="w3-container w3-third"),

            html.Div([
                dcc.Graph(
                    id="dotplot",
                    config = {"scrollZoom" : True, "displayModeBar" : True}
            )], id="dotplot-div", 
                style={"display" : "none"},
                className="w3-container w3-twothird w3-card"),
        ]),
        html.Div(style={"display" : "none"}, id="selected-read"),
        html.Div(style={"display" : "none"}, id="selected-ref"),
    ])

    @app.callback(
        Output("trackplot", "figure"),
        Output("info-table", "children"),
        Output("selection-card", "style"),
        Output("selected-ref", "children"),
        Output("selected-read", "children"),
        Input("trackplot-layer", "value"),
        Input("trackplot", "clickData"))
    def update_trackplot(layer, click):
        table = list()
        ref = aln = read = None
        card_style = {"display" : "none"}
        if click is not None:
            coord = click["points"][0]
            ref = coord["x"]

            if coord["curveNumber"] < len(tracks):
                track = tracks.all[coord["curveNumber"]]
                aln = track.alignments.iloc[coord["y"]]
                read = aln["read_id"]
                mref = tracks.coords.ref_to_mref(ref, aln.fwd)

                layers = track.layers.loc[(mref, aln.name)]["dtw"]

                table.append(html.Tr(html.Td(html.B("%s:%d" % (tracks.coords.ref_name, ref)), colSpan=2)))
                table.append(html.Tr(html.Td([html.B("Read "), read], colSpan=2)))
                for l in ["current", "dwell", "model_diff"]:
                    table.append(html.Tr([
                        html.Td(html.B(LAYERS["dtw"][l].label)), 
                        html.Td("%.3f"%layers[l], style={"text-align":"right"})]))

                card_style = {"display" : "block"}

        print(layer)

        layer = parse_layer(layer)

        print(layer)
        fig = Trackplot(
            tracks, layer, 
            select_ref=ref, select_read=read, 
            conf=conf).fig
        fig.update_layout(uirevision=True)

        print("DONE")

        return fig, table, card_style, ref, read

    @app.callback(
        #Output("dotplot", "figure"),
        Output("dotplot-div", "style"),
        State("selected-read", "children"),
        Input("dotplot-btn", "n_clicks"))
    def update_trackplot(read, n_clicks):
        if n_clicks is None:
            return {"display" : "hidden"}
        return {"display" : "block"}

    app.run_server(debug=True)
