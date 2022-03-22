"""Plot signal-to-reference alignment dotplots"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time

from .sigplot import Sigplot

from .. import config, nt
from ..dtw.aln_track import LAYERS, parse_layers
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..dtw.bcaln import Bcaln
from ..argparse import Opt, comma_split
from ..fast5 import parse_read_ids

class DotplotParams(config.ParamGroup):
    _name = "dotplot"
DotplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("bcaln_track", None, str, "Only display basecalled alignments from this track"),
    ("bcaln_error", False, bool, "Display basecalled alignment errors"),
    ("show_legend", True, bool, "Display legend"),
    ("select_ref", None, int, "Display a horizontal line at specified reference coordinate"),
    ("layers", [], None, ""),
)

class Dotplot:

    #TODO
    #plot_all(outfile) -> outfiles
    #iter_plots() -> figs
    #plot_read(read_id) -> fig

    REQ_LAYERS = [
        "start", "length", "middle", 
        "current", "dwell", "kmer", "base", 
        "bcaln.middle", "bcaln.error"
    ]

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("dotplot", *args, **kwargs)

        self.conf.tracks.layers = self.REQ_LAYERS + self.prms.layers

        if isinstance(self.prms.tracks, str) or self.prms.tracks is None:
            self.tracks = Tracks(conf=self.conf)
        elif isinstance(self.prms.tracks, Tracks):
            self.tracks = self.prms.tracks
        else:
            raise ValueError("Dotplot tracks parameter must be string or Tracks instance")

        self.layers = list(parse_layers(self.prms.layers, False))

        self.tracks.set_layers(self.REQ_LAYERS + self.layers)

        self.conf.load_config(self.tracks.conf)

        self.fig_config = {
                "toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None},
                "scrollZoom" : True, 
                "displayModeBar" : True}

    def iter_plots(self):
        t0 = time.time()
        for read_id, tracks in self.tracks.iter_reads():
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
            subplot_titles = [read_id], #+ None*(2*len(column_widths)-1)
            shared_xaxes=True,
            shared_yaxes=True)

        tracks_filter,colors_filter = zip(*[(t,c) for t,c in zip(tracks,self.conf.vis.track_colors) if t.name != self.prms.bcaln_track])
        #colors_filter = [t for t in tracks if t.name != self.prms.bcaln_track]

        Sigplot(tracks, conf=self.conf).plot(fig)

        hover_layers = [("dtw", "middle"),("dtw","kmer"),("dtw","current"),("dtw","dwell")] + self.layers
        #hover_layers += (l for l in self.prms.layers if l not in {"current","dwell"})
        hover_data = dict()

        coords = None

        for i,track in enumerate(tracks):

            track_hover = list()

            if track.coords is not None:
                coords = track.coords

            has_bcaln = "bcaln" in track.layers.columns.get_level_values("group")
            only_bcaln = self.prms.bcaln_track == track.name

            first_aln = True
            for aln_id, aln in track.alignments.iterrows():
                layers = track.get_aln_layers(aln_id)
                
                if has_bcaln:
                    self._plot_bcaln(fig, legend, layers)

                if not only_bcaln: 
                    fig.add_trace(go.Scattergl(
                        x=layers["dtw","start"], y=layers.index,
                        name=track.desc,
                        legendgroup=track.name,
                        line={"color":self.conf.vis.track_colors[i], "width":2, "shape" : "vh"},
                        hoverinfo="skip",
                        showlegend=first_aln
                    ), row=2, col=1)
                if only_bcaln: continue

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
                        #label = LAYERS["cmp"][stat].label
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

            hover_data[track.name] = pd.concat(track_hover)#.reset_index()

            if len(track_hover) > 0:
                hover_data[track.name] = track_hover[0]#.reset_index()

        if self.prms.bcaln_error:
            for i,track in enumerate(tracks):
                if not ("bcaln","error") in track.layers.columns: 
                    continue
                for aln_id, aln in track.alignments.iterrows():
                    self._plot_errors(fig, legend, track.get_aln_layers(aln_id))

        if len(hover_data) > 0:
            hover_data = pd.concat(hover_data, axis=1)
            hover_coords = hover_data.xs("middle", axis=1, level=2).mean(axis=1)

            hover_kmers = nt.kmer_to_str(
                hover_data.xs("kmer", 1, 2)
                          .fillna(method="pad", axis=1)
                          .iloc[:,-1])

            customdata = hover_data.drop(["kmer","middle"], axis=1, level=2).to_numpy()

            hover_rows = [
                "<b>" + coords.ref_name + ":%{y:,d} [%{text}]</b>"
            ]
            labels = [LAYERS[g][l].label for g,l in hover_layers[2:]]

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

        if not track.empty and track.all_fwd == self.conf.is_rna:
            fig.update_yaxes(autorange="reversed", row=2, col=1)
            fig.update_yaxes(autorange="reversed", row=2, col=2)

        if self.prms.select_ref is not None:
            fig.add_hline(y=self.prms.select_ref, line_color="red", row=2, col=1, opacity=0.5)
            i = hover_data.index.get_loc(self.prms.select_ref)

        strand = "+" if aln["fwd"] else "-"
        fig.update_yaxes(row=2, col=1,
            title_text=aln["ref_name"] + f" ({strand})")

        for i,(group,layer) in enumerate(self.layers):
            fig.update_xaxes(row=2, col=i+2,
                title_text=LAYERS[group][layer].label)

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
            margin={"t":50, "b":50},
            barmode="overlay",
            hoverdistance=20,
            dragmode="pan", 
            showlegend=self.prms.show_legend,
            legend={"bgcolor" : "#e6edf6"})#, scroll_zoom=True)

        return fig

    def _plot_bcaln(self, fig, legend, layers):
        fig.add_trace(go.Scattergl(
            x=layers["bcaln","start"], y=layers.index,#+2, #-1
            #x=layers["bcaln","start"], y=layers.index,#+2, #-1
            name="Basecalled Alignment",
            mode="markers", marker={"size":5,"color":"orange"},
            #line={"color":"orange", "width":2, "shape" : "vh"},
            legendgroup="bcaln",
            hoverinfo="skip",
            showlegend="bcaln_starts" not in legend
        ), row=2, col=1)
        legend.add("bcaln_starts")


    def _plot_errors(self, fig, legend, layers):
        if ("bcaln","error") not in layers.columns:
            return 

        errors = layers["bcaln","error"].dropna()
        sub = errors[errors.str.startswith("*")].str.slice(2)

        ins = errors[errors.str.startswith("+")]\
                    .str.slice(1)\
                    .map(list).explode()

        del_ = errors[errors.str.startswith("-")]\
               .str.slice(1)\
               .map(list).explode()

        #TODO global vis params
        #colors = ["#80ff80", "#6b93ff", "#ffbd00", "#ff8080"]
        linewidth = 3
        size = 15
        for b,base in enumerate(["a","c","g","t"]):
            refs = sub.index[sub.str.match(base)]
            if len(refs) > 0:
                fig.add_trace(go.Scattergl(
                    x=layers.loc[refs, ("bcaln","start")],
                    y=refs+2,
                    mode="markers",
                    marker_line_color=self.conf.vis.base_colors[b],
                    marker_line_width=linewidth,
                    marker_size=size,
                    marker_symbol="x-thin",
                    #hoverinfo="skip",
                    legendgroup="Bcaln Error",
                    name="SUB",
                    showlegend="bcaln_sub" not in legend,
                    legendrank=3
                ), row=2, col=1)
                legend.add("bcaln_sub")

            refs = ins.index[ins.str.match(base)]
            if len(refs) > 0:
                fig.add_trace(go.Scattergl(
                    x=layers.loc[refs, ("bcaln","start")],
                    y=refs+2,
                    mode="markers",
                    marker_line_color=self.conf.vis.base_colors[b],
                    marker_line_width=linewidth,
                    marker_size=size,
                    marker_symbol="cross-thin",
                    #hoverinfo="skip",
                    legendgroup="Bcaln Error",
                    name="INS",
                    showlegend="bcaln_ins" not in legend,
                    legendrank=4
                ), row=2, col=1)
                legend.add("bcaln_ins")

            refs = del_.index[del_.str.match(base)]
            if len(refs) > 0:
                starts = layers[("bcaln","start")].dropna()
                idxs = starts.index.get_indexer(refs,method="nearest")
                samps = starts.iloc[idxs]

                fig.add_trace(go.Scattergl(
                    x=samps,
                    y=refs+2,
                    mode="markers",
                    marker_line_color=self.conf.vis.base_colors[b],
                    marker_line_width=linewidth,
                    marker_size=size,
                    marker_symbol="line-ew",
                    #hoverinfo="skip",
                    name="DEL",
                    legendgroup="Bcaln Error",
                    showlegend="bcaln_del" not in legend,
                    legendrank=5
                ), row=2, col=1)
                legend.add("bcaln_del")


OPTS = (
    Opt("input", "tracks.io", nargs="+"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot.", required=True),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", "layers", type=comma_split),
    Opt(("-b", "--bcaln-track"), "dotplot"),
    Opt(("--multi-background"), "sigplot", action="store_true"),
    Opt(("--show-events"), "sigplot", action="store_true"),
    Opt(("--no-model"), "sigplot", action="store_true"),
    Opt(("--bcaln-error", "-e"), "dotplot", action="store_true"),
)

def main(conf):
    """Plot signal-to-reference alignment dotplots"""

    dotplots = Dotplot(conf=conf)
    for read_id, fig in dotplots.iter_plots():
        print(read_id)
        fig.write_html(
            conf.out_prefix + read_id + ".html", 
            config=dotplots.fig_config)
