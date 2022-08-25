import plotly.express as px
import pandas as pd
import numpy as np

import dash
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State
import sys

from .trackplot import Trackplot, PLOT_LAYERS
from .dotplot import Dotplot
from .refplot import Refplot
from .. import config
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..dtw.layers import LAYER_META, parse_layer
from ..argparse import Opt, comma_split
from ..config import Config
from ..fast5 import parse_read_ids


def browser(conf):
    """Interactive signal alignment genome browser"""
    conf.tracks.load_mat = True
    #conf.tracks.refstats_layers = ["cmp.mean_ref_dist"]
    conf.tracks.layers=["dtw","bcaln","cmp","bc_cmp"]
    sys.stderr.write("Loading tracks...\n")
    tracks = Tracks(conf=conf)
    sys.stderr.write("Starting server...\n")
    new_browser(tracks, conf)

def _icon_btn(icon, name=None, panel="", hide=False):
    style={"display" : "none" if hide else "inline-block"}
    if name is not None:
        href = f"javascript:{name}('{panel}')"
        id=f"{panel}-{name}"
    else:
        href = "javascript:void()"
        id=""

    return html.A(icon, id=id, className="material-icons w3-padding-24", href=href, style=style)

def _panel(title, name, content, settings=None, hide=False):
    style={"display" : "none" if hide else "block"}
    btns = list()

    if settings is not None:
        btns.append(_icon_btn("settings", "toggle_settings", name))
    btns += [
        _icon_btn("remove", "minimize", name),
        _icon_btn("add", "maximize", name, hide=True),
    ]

    ret = [html.Header(
        id=f"{name}-header", 
        className="w3-display-container w3-deep-purple", 
        style={"height":"40px"},
        children = [
            html.Div(
                html.H5(html.B(title)),
                className="w3-padding w3-display-left"),
            
            html.Div(children=btns, className="w3-display-right w3-padding"),
    ])]

    if settings is not None:
        ret.append(html.Div(
            settings, 
            id=f"{name}-settings", 
            style={"display" : "none"},
            className="w3-container w3-pale-blue"))
        
    ret.append(
        html.Div(content, id=f"{name}-body", className="w3-container"))

    return html.Div(
        html.Div(ret, id=f"{name}-card", className="w3-card"),
        id=f"{name}-panel", className="w3-panel", style=style)

def new_browser(tracks, conf):
    external_stylesheets = [
        "https://fonts.googleapis.com/icon?family=Material+Icons",
        "https://www.w3schools.com/w3css/4/w3.css",
    ]

    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
    app.title = "Uncalled4 Browser"

    layer_opts = [
        {"label" : LAYER_META.loc[(group,layer),"label"], "value" : f"{group}.{layer}"}
        for group,layer in tracks.aln_layers(PLOT_LAYERS)]
    
    track_opts = [
        {"label" : t.desc, "value" : t.name}
        for t in tracks.alns]

    app.layout = html.Div(children=[
        html.Div(
            html.H3(html.B("Uncalled4 Genome Browser")), 
            className="w3-container w3-deep-purple"),

        html.Div([
            html.Div(
                _panel("Trackplot", "trackplot", 
                    content=[
                        dcc.Dropdown(
                            options=layer_opts,
                            value=layer_opts[0]["value"], 
                            clearable=False, multi=False,
                            className="w3-padding",
                            id="trackplot-layer"),

                        dcc.Graph(#[dcc.Loading(type="circle"),
                            id="trackplot",
                            config = {"scrollZoom" : True, "displayModeBar" : True})

                    ], settings=[
                        dcc.Dropdown(
                            options=track_opts,
                            value=tracks.track_names,
                            clearable=False, multi=True, searchable=False,
                            className="w3-padding",
                            id="track-dropdown"),
                        dcc.Checklist(
                            id="trackplot-checklist",
                            className="w3-container w3-padding",
                            labelStyle={"display" : "block"},
                            inputClassName="w3-padding",
                            options=[
                                {"label" : "Show legend", "value" : "show_legend"},
                                {"label" : "Full overlap", "value" : "full_overlap"},
                                {"label" : "Shared reads only", "value" : "share_reads"},
                            ], value=["show_legend"])
                ]) #end trackplot panel
            , className="w3-half"), #end trackplot div

            html.Div([
                _panel("Selection", "selection",
                        html.Table([], id="info-table")),

                _panel("Refplot", "refplot",
                    content=dcc.Graph(
                        id="refplot",
                        config = {"scrollZoom" : True, "displayModeBar" : True}
                    ), settings=[
                        dcc.Checklist(
                            id="refplot-checklist",
                            className="w3-container w3-padding",
                            labelStyle={"display" : "block"},
                            inputClassName="w3-padding",
                            options=[
                                {"label" : "Show legend", "value" : "show_legend"},
                                {"label" : "Show model current", "value" : "show_model"},
                                {"label" : "Always color bases", "value" : "multi_background"},
                            ], value=["show_legend", "show_model"])
                    ], hide=True),

                _panel("Dotplot", "dotplot",
                    content=dcc.Graph(
                        id="dotplot",
                        config = {"scrollZoom" : True, "displayModeBar" : True}
                    ), settings=[
                        dcc.Checklist(
                            id="dotplot-checklist",
                            className="w3-container w3-padding",
                            labelStyle={"display" : "block"},
                            inputClassName="w3-padding",
                            options=[
                                {"label" : "Show legend", "value" : "show_legend"},
                                {"label" : "Show model current", "value" : "show_model"},
                                {"label" : "Always color bases", "value" : "multi_background"},
                            ], value=["show_legend", "show_model"])
                    ], hide=True),
            ], className="w3-half"),
            html.Div([
            ], className="w3-full"),

        ]),
        html.Div(style={"display" : "none"}, id="selected-read"),
        html.Div(style={"display" : "none"}, id="selected-ref"),
        html.Div(style={"display" : "none"}, id="read-changed"),
    ])

    @app.callback(
        Output("trackplot", "figure"),
        Output("info-table", "children"),
        Output("selection-panel", "style"),
        Output("selected-ref", "children"),
        Output("selected-read", "children"),
        Output("read-changed", "children"),
        Input("track-dropdown", "value"),
        Input("trackplot-checklist", "value"),
        Input("trackplot-layer", "value"),
        Input("trackplot", "clickData"),
        State("selected-read", "children"))
    def update_trackplot(track_names, checklist, layer, click, prev_read):
        table = list()
        ref = aln = read = None
        card_style = {"display" : "none"}

        checklist = set(checklist)

        chunk = tracks
        uirev = True

        shared_reads = "share_reads" in checklist
        full_overlap = "full_overlap" in checklist

        if len(track_names) == 0:
            track_names = None

        chunk = chunk.slice(tracks=track_names, shared_reads=shared_reads, full_overlap=full_overlap)

        if click is not None:
            coord = click["points"][0]
            ref = coord["x"]

            if coord["curveNumber"] < len(chunk):
                track = chunk.alns[coord["curveNumber"]]
                aln = track.alignments.iloc[coord["y"]]
                read = aln["read_id"]

                layers = track.layers.loc[(ref, aln.name)]["dtw"]

                table.append(html.Tr(html.Td(html.B("%s:%d" % (chunk.coords.ref_name, ref)), colSpan=2)))
                table.append(html.Tr(html.Td([html.B("Read "), read], colSpan=2)))
                for l in ["current", "dwell", "model_diff", "events"]:
                    table.append(html.Tr([
                        html.Td(html.B(LAYER_META.loc[("dtw",l), "label"])), 
                        html.Td("%.3f"%layers[l], style={"text-align":"right"})]))

                card_style = {"display" : "block"}

        layer, = parse_layer(layer)

        fig = Trackplot(
            chunk, [("mat", layer)], 
            select_ref=ref, select_read=read, 
            share_reads="share_reads" in checklist,
            show_legend="show_legend" in checklist,
            conf=conf).fig
        fig.update_layout(uirevision=uirev)

        return fig, table, card_style, ref, read, (read != prev_read)

    @app.callback(
        Output("dotplot", "figure"),
        Output("dotplot-panel", "style"),
        #State("selected-read", "children"),
        #Input("dotplot-btn", "n_clicks"))
        Input("track-dropdown", "value"),
        Input("dotplot-checklist", "value"),
        Input("trackplot-layer", "value"),
        Input("selected-ref", "children"),
        Input("selected-read", "children"),
        Input("read-changed", "children"))
    def update_dotplot(track_names, flags, layer, ref, read, read_changed):
        if read is None:
            return {}, {"display" : "hidden"}

        flags = set(flags)

        conf = Config(conf=tracks.conf)
        conf.sigplot.multi_background="multi_background" in flags
        conf.sigplot.no_model="show_model" not in flags

        chunk = tracks.slice(tracks=track_names if len(track_names) > 0 else None)

        fig = Dotplot(
            chunk, 
            select_ref=ref, 
            show_legend="show_legend" in flags,
            layers=list(parse_layer(layer)),
            conf=conf).plot(read)
        #fig.update_layout(uirevision=read_changed)

        return fig, {"display" : "block"}


    @app.callback(
        Output("refplot", "figure"),
        Output("refplot-panel", "style"),
        #State("selected-read", "children"),
        #Input("dotplot-btn", "n_clicks"))
        Input("track-dropdown", "value"),
        Input("dotplot-checklist", "value"),
        Input("trackplot-layer", "value"),
        Input("selected-ref", "children"),
        Input("selected-read", "children"),
        Input("read-changed", "children"))
    def update_refplot(track_names, flags, layer, ref, read, read_changed):
        if read is None:
            return {}, {"display" : "hidden"}

        flags = set(flags)

        conf = Config(conf=tracks.conf)
        conf.sigplot.multi_background="multi_background" in flags
        conf.sigplot.no_model="show_model" not in flags

        chunk = tracks.slice(tracks=track_names if len(track_names) > 0 else None)

        fig = Refplot(
            chunk, 
            layer=layer,
            kmer_coord=ref,
            conf=conf).fig
        #fig.update_layout(uirevision=read_changed)

        return fig, {"display" : "block"}

    app.run_server(port=conf.browser_port, debug=True)
