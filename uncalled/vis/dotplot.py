"""Plot signal-to-reference alignment dotplots"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from .sigplot import Sigplot

from .. import config, nt
from ..dtw.aln_track import LAYERS, parse_layers
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..dtw.bcaln import Bcaln
from ..argparse import Opt, comma_split
from ..fast5 import parse_read_ids
from ..sigproc import ProcRead

track_colors = ["#AA0DFE", "#1CA71C", "#6A76FC"]

class DotplotParams(config.ParamGroup):
    _name = "dotplot"
DotplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("layers", [], None, ""),
    ("track_colors", ["#AA0DFE", "#1CA71C", "#4676FF"], list, ""),
)

class Dotplot:

    #TODO
    #plot_all(outfile) -> outfiles
    #iter_plots() -> figs
    #plot_read(read_id) -> fig

    _req_layers = [
        "start", "length", "middle", 
        "current", "dwell", "kmer", "base", 
        "bcaln.start"#, "cmp.mean_ref_dist"
    ]

    def __init__(self, *args, **kwargs):
        conf, self.prms = config._init_group("dotplot", *args, **kwargs)

        if isinstance(self.prms.tracks, str) or self.prms.tracks is None:
            self.tracks = Tracks(conf=conf)
        elif isinstance(self.prms.tracks, Tracks):
            self.tracks = self.prms.tracks
        else:
            raise ValueError("Dotplot tracks parameter must be string or Tracks instance")

        self.conf = self.tracks.conf

        self.prms.layers = list(parse_layers(self.prms.layers, False))

        self.tracks.set_layers(self._req_layers + conf.dotplot.layers)

        self.fig_config = {
                "toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None},
                "scrollZoom" : True, 
                "displayModeBar" : True}

    def iter_plots(self):
        for read_id, tracks in self.tracks.iter_reads():
            print(read_id)
            yield read_id, self._plot(read_id, tracks)

    def plot(self, read_id):
        for read_id, tracks in self.tracks.iter_reads([read_id]):
            return self._plot(read_id, tracks)

    def show(self, read_id):
        self.plot(read_id).show(config=self.fig_config)

    def _plot(self, read_id, tracks):

        #if tracks[0].has_group("cmp"):
        #    cmp_stats = list(tracks[0].layers["cmp"].columns)
        #else:
        cmp_stats = []

        column_widths=[6]+[1]*(len(self.prms.layers)+len(cmp_stats))

        fig = make_subplots(
            rows=2, cols=len(column_widths), 
            row_heights=[1,3],
            column_widths=column_widths,
            vertical_spacing=0.01,
            horizontal_spacing=0.01,
            subplot_titles = [read_id], #+ None*(2*len(column_widths)-1)
            shared_xaxes=True,
            shared_yaxes=True)

        Sigplot(tracks, track_colors=self.prms.track_colors, conf=self.conf).plot(fig)

        hover_layers = [("dtw", "middle"),("dtw","kmer"),("dtw","current"),("dtw","dwell")] + self.prms.layers
        #hover_layers += (l for l in self.prms.layers if l not in {"current","dwell"})
        hover_data = dict()

        for i,track in enumerate(tracks):

            track_hover = list()

            has_bcaln = "bcaln" in track.layers.columns.get_level_values("group")

            first_aln = True
            for aln_id, aln in track.alignments.iterrows():
                layers = track.get_aln_layers(aln_id)

                track_hover.append(layers[hover_layers])
                
                if has_bcaln:
                    fig.add_trace(go.Scatter(
                        x=layers["bcaln","start"], y=layers.index+Bcaln.K-2,
                        name="Basecalled Alignment",
                        mode="markers", marker={"size":5,"color":"orange"},
                        legendgroup="bcaln",
                        hoverinfo="skip",
                        showlegend=first_aln
                    ), row=2, col=1)

                fig.add_trace(go.Scatter(
                    x=layers["dtw","start"], y=layers.index,
                    name=track.desc,
                    legendgroup=track.name,
                    line={"color":self.prms.track_colors[i], "width":2, "shape" : "hv"},
                    hoverinfo="skip",
                    showlegend=first_aln
                ), row=2, col=1)
                    

                first_aln = False

                for j,layer in enumerate(self.prms.layers):
                    if layer[0] != "cmp":
                        fig.add_trace(go.Scatter(
                            x=layers[layer], y=layers.index-0.5,
                            name=track.desc, 
                            line={
                                "color" : self.prms.track_colors[i], 
                                "width":2, "shape" : "hv"},
                            legendgroup=track.name, showlegend=False,
                        ), row=2, col=j+2)

                    elif len(layers[layer].dropna()) > 0:
                        color= "rgba(255,0,0,1)"
                        #label = LAYERS["cmp"][stat].label
                        fig.add_trace(go.Bar(
                            x=track.layers[layer], 
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

            #hover_data[track.name] = pd.concat(track_hover)#.reset_index()
            hover_data[track.name] = track_hover[0]#.reset_index()

        hover_data = pd.concat(hover_data, axis=1)
        hover_coords = hover_data.xs("middle", axis=1, level=2).mean(axis=1)

        hover_kmers = nt.kmer_to_str(
            hover_data.xs("kmer", 1, 2)
                      .fillna(method="pad", axis=1)
                      .iloc[:,-1])


        customdata = hover_data.drop(["kmer","middle"], axis=1, level=2).to_numpy()

        hover_rows = [
            "<b>" + track.coords.ref_name + ":%{y:,d} [%{text}]</b>"
        ]
        labels = [LAYERS[g][l].label for g,l in hover_layers[2:]]

        for i,label in enumerate(labels):
            s = label
            fields = list()
            for j in range(len(tracks)):
                k = len(labels) * j + i
                fields.append(
                    '<span style="color:%s;float:right"><b>%%{customdata[%d]:.2f}</b></span>' % 
                    (self.prms.track_colors[j], k))
            hover_rows.append(s + ": " + ", ".join(fields))


        fig.add_trace(go.Scatter(
            x=hover_coords, y=hover_data.index,
            mode="markers", marker={"size":0,"color":"rgba(0,0,0,0)"},
            name="",
            customdata=customdata,
            hovertemplate="<br>".join(hover_rows),
            hoverlabel={"bgcolor":"rgba(255,255,255,1)"},
            text=hover_kmers,
            showlegend=False
        ), row=2,col=1)
            #customdata=hoverdata.to_numpy(),

        if track.coords.fwd == self.conf.is_rna:
            fig.update_yaxes(autorange="reversed", row=2, col=1)
            fig.update_yaxes(autorange="reversed", row=2, col=2)

        fig.update_yaxes(row=2, col=1,
            title_text="Reference (%s)" % aln["ref_name"])

        for i,(group,layer) in enumerate(self.prms.layers):
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
            barmode="overlay",
            hoverdistance=20,
            dragmode="pan", 
            legend={"bgcolor" : "#e6edf6"})#, scroll_zoom=True)

        return fig


OPTS = (
    Opt("input", "tracks", nargs="+"),
    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot.", required=True),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", "layers", type=comma_split),
)

def main(conf):
    """Plot signal-to-reference alignment dotplots"""

    dotplots = Dotplot(conf=conf)
    for read_id, fig in dotplots.iter_plots():
        fig.write_html(
            conf.out_prefix + read_id + ".html", 
            config=dotplots.fig_config)
