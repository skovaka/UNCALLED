"""Plot signal-to-reference alignment dotplots"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time
import sys

from .sigplot import Sigplot

from .. import config
from ..aln import LAYER_META, parse_layers
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..argparse import Opt, comma_split

class Dotplot:

    #TODO
    #plot_all(outfile) -> outfiles
    #iter_plots() -> figs
    #plot_read(read_id) -> fig

    REQ_LAYERS = [
        "dtw.start", "dtw.length", "dtw.middle_sec", 
        "dtw.current", "dtw.current_sd", "dtw.dwell", "seq.kmer",  
        "moves.middle_sec", "seq.current",
        "dtw.start_sec", "dtw.length_sec"
    ]

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("dotplot", *args, **kwargs)

        req_layers = self.REQ_LAYERS.copy()
        #if self.prms.show_bands:
        #    req_layers.append("band")

        self.conf.tracks.layers = req_layers + self.prms.layers

        if isinstance(self.prms.tracks, str) or self.prms.tracks is None:
            self.tracks = Tracks(conf=self.conf)
        elif isinstance(self.prms.tracks, Tracks):
            self.tracks = self.prms.tracks
        else:
            raise ValueError("Dotplot tracks parameter must be string or Tracks instance")

        self.layers = list(parse_layers(self.prms.layers))

        self.tracks.set_layers(req_layers + self.layers)

        self.conf.load_config(self.tracks.conf)

        self.fig_config = {
                "toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None},
                "scrollZoom" : True, 
                "displayModeBar" : True}

    def iter_plots(self):
        t0 = time.time()
        #for read_id, aln in self.tracks.iter_reads(return_tracks=True):
        for read_id, tracks in self.tracks.iter_reads(return_tracks=True):
            #print([a for a in aln])
            #print(aln.to_pandas(["dtw"]))
            #print(read_id, aln)
            yield read_id, self._plot(read_id, tracks)
            t0 = time.time()

    def plot(self, read_id):
        chunk = self.tracks.slice(reads=[read_id])
        return self._plot(read_id, chunk)
        #for read_id, tracks in self.tracks.iter_reads_([read_id]):
        #    return self._plot(read_id, tracks)

    def show(self, read_id):
        self.plot(read_id).show(config=self.fig_config)

    def _plot(self, read_id, tracks):
        #if tracks[0].has_group("cmp"):
        #    cmp_stats = list(tracks[0].layers["cmp"].columns)
        #else:
        cmp_stats = []

        column_widths=[6]+[1]*(len(self.layers)+len(cmp_stats))

        legend = set()

        fig = make_subplots(
            rows=2, cols=len(column_widths), 
            row_heights=[1,3],
            column_widths=column_widths,
            vertical_spacing=0.01,
            horizontal_spacing=0.01,
            #subplot_titles = [read_id], #+ None*(2*len(column_widths)-1)
            shared_xaxes=True,
            shared_yaxes=True)

        tracks_filter,colors_filter = zip(*[(t,c) for t,c in zip(tracks.alns,self.conf.vis.track_colors) if t.name != self.prms.moves_track])
        #colors_filter = [t for t in tracks if t.name != self.prms.moves_track]

        Sigplot(tracks, conf=self.conf).plot(fig)

        hover_layers = [("dtw", "middle_sec"),("seq","kmer"),("dtw","current"),("dtw","dwell")] + self.layers
        #hover_layers += (l for l in self.prms.layers if l not in {"current","dwell"})
        hover_data = dict()

        coords = tracks.coords

        flipped = True
        fwd = True
        moves_plotted = False

        active_tracks = tracks.alignments.index.unique("track")
        print(tracks.alignments)

        for i,track in enumerate(tracks.alns):
            print(track.name)
            if track.name not in active_tracks: continue

            track_hover = list()

            has_moves = "moves" in tracks.layers.columns.get_level_values("group")
            only_moves = self.prms.moves_track == track.name
            model = track.model

            first_aln = True
            for aln_id, aln in tracks.alignments.loc[track.name].iterrows():
                layers = tracks.layers \
                              .loc[(track.name,aln_id)] \
                              .reset_index() \
                              .set_index("seq.pos") \
                              .sort_values(("dtw","start_sec"))
                              #droplevel("seq.fwd") #,slice(None)), slice(None)] \
                              #.droplevel("aln.id")
                
                fwd = fwd and aln["fwd"]
                flipped = flipped and aln["fwd"] == self.conf.is_rna

                if self.prms.show_bands and "band" in layers:
                    bands = layers["band"].dropna(how="any")
                    fig.add_trace(go.Scattergl(
                        x=bands["sample_start"], 
                        y=bands.index,
                        line={"color" : "orange"},
                        fillcolor="orange",
                        opacity=0.2,
                        fill="tonexty",
                        name="DTW Band",
                        legendgroup="band",
                    ), row=2, col=1)
                    fig.add_trace(go.Scattergl(
                        x=bands["sample_end"], 
                        y=bands["ref_end"],
                        line={"color" : "orange"},
                        fillcolor="orange",
                        opacity=0.2,
                        fill="tonexty",
                        legendgroup="band",
                        showlegend=False
                    ), row=2, col=1)
                
                if has_moves and not moves_plotted:
                    moves_plotted = True
                    self._plot_moves(fig, legend, layers)

                if not only_moves: 
                    fig.add_trace(go.Scattergl(
                        x=layers["dtw","start_sec"], y=layers.index,
                        name=track.desc,
                        legendgroup=track.name,
                        line={
                            "color":self.conf.vis.track_colors[i], 
                            "width":2, "shape" : "hv" }, #if not flipped else "vh" },
                        hoverinfo="skip",
                        showlegend=first_aln
                    ), row=2, col=1)
                if only_moves: continue

                track_hover.append(layers[hover_layers])
                    

                first_aln = False

                for j,layer in enumerate(self.layers):
                    if layer[0] != "cmp":
                        fig.add_trace(go.Scattergl(
                            x=layers[layer], y=layers.index+0.5,
                            name=track.desc, 
                            line={
                                "color" : self.conf.vis.track_colors[i], 
                                "width":2, "shape" : "hv"},
                            legendgroup=track.name, showlegend=False,
                        ), row=2, col=j+2)

                    elif layer in layers and len(layers[layer].dropna()) > 0:
                        color= "rgba(255,0,0,1)"
                        fig.add_trace(go.Bar(
                            x=track.layers[layer].fillna(1.0), 
                            y=track.layer_refs, #TODO try vhv
                            base=0,
                            name="DTW Compare",
                            orientation="h",
                            width=1,
                            marker={"color":color,"line":{"color":color,"width":0.5}},
                            legendgroup="cmp",
                            #showlegend=i==0
                        ), row=2, col=j+2)
            #fig.update_xaxes(row=2, col=cmp_col + i,
            #    title_text=label)


            if len(track_hover) > 0:
                hover_data[track.name] = pd.concat(track_hover)#.reset_index()
                hover_data[track.name] = track_hover[0]#.reset_index()

        #if self.prms.moves_error:
        #    for i,track in enumerate(tracks):
        #        if not ("moves","error") in track.layers.columns: 
        #            continue
        #        for aln_id, aln in track.alignments.iterrows():
        #            layers = track.layers.loc[(slice(None),aln_id)].droplevel("aln.id")
        #            self._plot_errors(fig, legend, layers)

        if len(hover_data) > 0:
            hover_data = pd.concat(hover_data, axis=1)
            hover_coords = hover_data.xs("middle_sec", axis=1, level=2).mean(axis=1)

            hover_kmers = model.kmer_to_str(
                hover_data.xs("kmer", 1, 2)
                          .fillna(method="pad", axis=1)
                          .iloc[:,-1])

            customdata = hover_data.drop(["kmer","middle_sec"], axis=1, level=2).to_numpy()

            hover_rows = [
                "<b>" + coords.name + ":%{y:,d} [%{text}]</b>"
            ]
            labels = [LAYER_META.loc[(g,l),"label"] for g,l in hover_layers[2:]]

            for i,label in enumerate(labels):
                s = label
                fields = list()
                for j in range(len(tracks_filter)):
                    k = len(labels) * j + i
                    fields.append(
                        '<span style="color:%s;float:right"><b>%%{customdata[%d]:.2f}</b></span>' % 
                        (self.conf.vis.track_colors[j], k))
                hover_rows.append(s + ": " + ", ".join(fields))


            fig.add_trace(go.Scattergl(
                x=hover_coords, y=hover_data.index,
                mode="markers", marker={"size":0,"color":"rgba(0,0,0,0)"},
                name="",
                customdata=customdata,
                hovertemplate="<br>".join(hover_rows),
                hoverlabel={"bgcolor":"rgba(255,255,255,1)"},
                text=hover_kmers,
                showlegend=False
            ), row=2,col=1)

        #if not track.empty and track.all_fwd == self.conf.is_rna:
        if flipped:
            fig.update_yaxes(autorange="reversed", row=2, col=1)
            fig.update_yaxes(autorange="reversed", row=2, col=2)

        if self.prms.select_ref is not None:
            fig.add_hline(y=self.prms.select_ref, line_color="red", row=2, col=1, opacity=0.5)
            i = hover_data.index.get_loc(self.prms.select_ref)

        strand = "+" if fwd else "-"
        fig.update_yaxes(row=2, col=1,
            title_text=tracks.coords.name + f" ({strand})")

        for i,(group,layer) in enumerate(self.layers):
            fig.update_xaxes(row=2, col=i+2,
                title_text=LAYER_META.loc[(group,layer),"label"])

        axis_kw = dict(
            showspikes=True,
            spikemode="across",
            spikecolor="darkgray",
            spikethickness=1)

        fig.update_xaxes(**axis_kw)
        fig.update_yaxes(**axis_kw)
        fig.update_yaxes(
            tickformat=",d", tickangle=-90, nticks=3
        )#showticklabels=False)
        #fig.update_yaxes(showspikes=True)

        fig.update_xaxes(
            title_text="Time (s)", 
            #tickformat=".3f", 
            #showspikes=True,
            row=2, col=1)

        fig.update_layout(
            #hovermode="x unified",
            margin={"l":100,"r":50},#, "b":50},
            barmode="overlay",
            hoverdistance=20,
            dragmode="pan", 
            showlegend=self.prms.show_legend,
            legend={
                "y" : 1.05, "yanchor" : "bottom",
                #"x" : 0, "xanchor" : "left",
                "traceorder" : "normal",
                "orientation" : "h", "bgcolor" : "#e6edf6"})#, scroll_zoom=True)

        return fig

    def _plot_moves(self, fig, legend, layers):
        fig.add_trace(go.Scattergl(
            x=layers["moves","middle_sec"], y=layers.index,#-2, #-1
            #x=layers["moves","start_sec"], y=layers.index,#+2, #-1
            name="Basecalled Alignment",
            mode="markers", marker={"size":5,"color":"orange"},
            #line={"color":"orange", "width":2, "shape" : "vh"},
            legendgroup="moves",
            hoverinfo="skip",
            showlegend="moves_starts" not in legend
        ), row=2, col=1)
        legend.add("moves_starts")


    #def _plot_errors(self, fig, legend, layers):
    #    if ("moves","error") not in layers.columns:
    #        return 

    #    errors = layers["moves","error"].dropna()
    #    sub = errors[errors.str.startswith("*")].str.slice(2)

    #    ins = errors[errors.str.startswith("+")]\
    #                .str.slice(1)\
    #                .map(list).explode()

    #    del_ = errors[errors.str.startswith("-")]\
    #           .str.slice(1)\
    #           .map(list).explode()

    #    #TODO global vis params
    #    #colors = ["#80ff80", "#6b93ff", "#ffbd00", "#ff8080"]
    #    linewidth = 3
    #    size = 15
    #    for b,base in enumerate(["a","c","g","t"]):
    #        refs = sub.index[sub.str.match(base)]
    #        if len(refs) > 0:
    #            fig.add_trace(go.Scattergl(
    #                x=layers.loc[refs, ("moves","start")],
    #                y=refs+2,
    #                mode="markers",
    #                marker_line_color=self.conf.vis.base_colors[b],
    #                marker_line_width=linewidth,
    #                marker_size=size,
    #                marker_symbol="x-thin",
    #                #hoverinfo="skip",
    #                legendgroup="Bcaln Error",
    #                name="SUB",
    #                showlegend="moves_sub" not in legend,
    #                legendrank=3
    #            ), row=2, col=1)
    #            legend.add("moves_sub")

    #        refs = ins.index[ins.str.match(base)]
    #        if len(refs) > 0:
    #            fig.add_trace(go.Scattergl(
    #                x=layers.loc[refs, ("moves","start")],
    #                y=refs+2,
    #                mode="markers",
    #                marker_line_color=self.conf.vis.base_colors[b],
    #                marker_line_width=linewidth,
    #                marker_size=size,
    #                marker_symbol="cross-thin",
    #                #hoverinfo="skip",
    #                legendgroup="Bcaln Error",
    #                name="INS",
    #                showlegend="moves_ins" not in legend,
    #                legendrank=4
    #            ), row=2, col=1)
    #            legend.add("moves_ins")

    #        refs = del_.index[del_.str.match(base)]
    #        if len(refs) > 0:
    #            starts = layers[("moves","start")].dropna()
    #            idxs = starts.index.get_indexer(refs,method="nearest")
    #            samps = starts.iloc[idxs]

    #            fig.add_trace(go.Scattergl(
    #                x=samps,
    #                y=refs+2,
    #                mode="markers",
    #                marker_line_color=self.conf.vis.base_colors[b],
    #                marker_line_width=linewidth,
    #                marker_size=size,
    #                marker_symbol="line-ew",
    #                #hoverinfo="skip",
    #                name="DEL",
    #                legendgroup="Bcaln Error",
    #                showlegend="moves_del" not in legend,
    #                legendrank=5
    #            ), row=2, col=1)
    #            legend.add("moves_del")



def dotplot(conf):
    """Plot signal-to-reference alignment dotplots"""

    conf.tracks.load_fast5 = True
    dotplots = Dotplot(conf=conf)
    save = conf.out_prefix is not None
    for read_id, fig in dotplots.iter_plots():
        if save:
            fig.write_html(
                conf.out_prefix + read_id + ".html", 
                config=dotplots.fig_config)
        else:
            fig.show(config=dotplots.fig_config)
            sys.stderr.flush()
            sys.stdout.write("Press enter to plot next read, or type \"exit\"\n")
            sys.stdout.flush()
            choice = sys.stdin.readline().strip().lower()
            if choice == "exit": break
