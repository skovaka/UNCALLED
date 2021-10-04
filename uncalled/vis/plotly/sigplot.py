import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from ... import nt

from ... import config
from ...dtw.track import LAYER_META
from ...index import str_to_coord
from ...dtw.track_io import TrackIO
from ...argparse import Opt, comma_split
from ...fast5 import parse_read_ids
from ...sigproc import ProcRead

class SigplotParams(config.ParamGroup):
    _name = "sigplot"
SigplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    #("ref_bounds", None, str_to_coord, "DTW aligment tracks"),
    #("reads", None, None, "Reads to plot"),
    ("track_colors", ["purple", "darkgreen", "royalblue", "crimson"], list, ""),
    ("fill_layer", "base", str, ""),
    ("fill_track", 0, None, "")
)

class Sigplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("sigplot", *args, **kwargs)
        self.tracks = self.prms.tracks #TODO do this better

        alns = self.tracks[0].alignments.index
        for t in self.tracks:
            alns = alns.union(t.alignments.index)
        self.alns = alns.sort_values()

    def plot(self, fig=None, row=1, col=1):
        if fig is None:
            fig = make_subplots(
                rows=len(self.alns), cols=1, 
                vertical_spacing=0.01,
                shared_xaxes=True)

        for aln_id in self.alns:
            self.plot_aln(fig, aln_id, row, col)
            row += 1

        fig.update_layout(yaxis={"fixedrange" : True}, dragmode="pan")

    def plot_aln(self, fig, aln_id, row=1, col=1):
        for i,track in enumerate(self.tracks):
            track_color = self.prms.track_colors[i]
            
            if aln_id not in track.alignments.index:
                print(aln_id)
                print(track.alignments.index)
                continue

            aln = track.alignments.loc[aln_id]

            read = ProcRead(track.fast5s[aln["read_id"]], conf=self.conf)

            dtw = track.layers["dtw"].loc[aln_id]

            samp_min = dtw["start"].min()
            max_i = dtw["start"].argmax()
            samp_max = dtw["start"].iloc[max_i] + dtw["length"].iloc[max_i]
            samps = np.arange(samp_min, samp_max)
            raw_norm = read.get_norm_signal(samp_min, samp_max)
            mask = ((raw_norm >= 40) &
                    (raw_norm <= self.conf.event_detector.max_mean))

            kmers = track.coords.kmers[dtw.index]
            model_current = track.model[kmers]

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
                ), row=row, col=col)

            fig.add_trace(go.Scattergl(
                x=samps[mask], y=raw_norm[mask],
                hoverinfo="skip",
                name="Raw Signal",
                mode="markers",
                marker={"size":2, "color":"black"}
            ), row=row, col=col)

            fig.add_trace(go.Scattergl(
                name = "Model Current",
                x=dtw["start"], y=model_current,
                line={"color":"white", "width":2, "shape" : "hv"}
            ), row=row, col=col)


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
        fast5_read = fast5s[read_id]
        if isinstance(fast5_read, ProcRead):
            read = fast5_read
        else:
            read = ProcRead(fast5_read, conf=io.conf)

        fig = dotplot(conf, read, tracks)

        #fig.show()

        fig.write_html(conf.out_prefix + read_id + ".html", config={"scrollZoom" : True, "displayModeBar" : True})
