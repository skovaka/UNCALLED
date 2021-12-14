import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from .. import nt

from .. import config
from ..dtw.aln_track import LAYERS
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..argparse import Opt, comma_split
from ..fast5 import parse_read_ids

class SigplotParams(config.ParamGroup):
    _name = "sigplot"
SigplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    #("ref_bounds", None, str_to_coord, "DTW aligment tracks"),
    #("reads", None, None, "Reads to plot"),
    ("max_reads", 10, int, ""),
    ("no_model", False, bool, "Will not plot the expected reference signal if True"),
    ("multi_background", False, bool, "Will plot multiple stacked background colors for multiple tracks if True"),
    ("yaxis_fixed", True, bool, ""),
    ("track_colors", ["purple", "darkgreen", "royalblue", "crimson"], list, ""),
    #("base_colors", ["#80ff80", "#6b93ff", "#ffbd00", "#ff8080"], list, "Colors for each base (A,C,G,T/U)"), 
    ("fill_layer", "base", str, ""),
    ("fill_track", 0, None, "")
)

class Sigplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("sigplot", *args, **kwargs)
        self.tracks = self.prms.tracks #TODO do this better

        self._legend = set()

        #reads = pd.Index([])
        #for t in self.tracks:
        #    reads = reads.union(t.alignments["read_id"]).unique()
        #if len(reads) > self.prms.max_reads:
        #    reads = np.random.choice(reads, self.prms.max_reads, False)
        self.reads = np.sort(self.tracks.get_all_reads())
        
    def plot(self, fig=None, row=1, col=1):
        if fig is None:
            fig = make_subplots(
                rows=len(self.reads), cols=1, 
                vertical_spacing=0.01)
            fig.update_yaxes(fixedrange=True)
            fig.update_layout(dragmode="pan")

        for read_id in self.reads:
            self.plot_read(fig, read_id, row, col)
            row += 1
        
        return fig

            
    def _plot_bases(self, fig, dtw, ymin, ymax, row, col):
        bases = nt.kmer_base(dtw["kmers"], 2)
        for base, color in enumerate(self.conf.vis.base_colors):
            base_dtw = dtw[bases == base]
            starts = base_dtw['start']
            ends = starts + base_dtw['length'] - 1
            nones = [None]*len(base_dtw)

            ys = [ymax,ymax,ymin,ymin,None]*len(base_dtw)
            xs = np.dstack([starts, ends, ends, starts, nones]).reshape(len(ys))

            fig.add_trace(go.Scattergl(
                x=xs,y=ys, fill="toself",
                fillcolor=color,
                hoverinfo="skip",
                mode="none",
                legendgroup="bases",
                showlegend= color not in self._legend,
                line={"width" : 0},
                name=nt.base_to_char(base),
                legendrank=1
            ), row=row, col=col)
            self._legend.add(color)

    def plot_read(self, fig, read_id, row=1, col=1):

        track_dtws  = list()
        
        current_min = samp_min = np.inf
        current_max = samp_max = 0
        for i,track in enumerate(self.tracks.alns):
            track_color = self.prms.track_colors[i]

            alns = track.alignments.query("@read_id == read_id")
            aln_ids = alns.index

            dtws = list()
            
            for aln_id,aln in alns.iterrows():
                dtw = track.layers["dtw"].xs(aln_id, level="aln_id")
                #dtw["kmers"] = track.coords.ref_kmers.loc[(aln.fwd, dtw.index)]
                dtw["kmers"] = track.kmers.xs(aln_id, level="aln_id")
                dtw["model_current"] = track.model[dtw["kmers"]]
                dtws.append(dtw)

                max_i = dtw["start"].argmax()
                samp_min = min(samp_min, dtw["start"].min())
                samp_max = max(samp_max, dtw["start"].iloc[max_i] + dtw["length"].iloc[max_i])

                current_min = min(current_min, dtw["model_current"].min())
                current_max = max(current_max, dtw["model_current"].max())

            if len(dtws) > 0:
                track_dtws.append(pd.concat(dtws).sort_index())

        read = ProcRead(track.fast5s[read_id], conf=self.conf)
        signal = read.get_norm_signal(samp_min, samp_max)

        #TODO set global signal min/max
        mask = ((signal >= 40) &
                (signal <= self.conf.event_detector.max_mean))
        
        samples = np.arange(samp_min, samp_max)[mask]
        signal = signal[mask]

        current_min = min(current_min, signal.min())
        current_max = max(current_max, signal.max())

        if len(self.tracks) == 1:
            self._plot_bases(fig, track_dtws[0], current_min, current_max, row, col)
            colors = ["white"]
            dtw_kws = [{}]
        else:

            if self.prms.multi_background:
                ys = np.linspace(current_min, current_max, len(self.tracks)+1)
                dy = (ys[1]-ys[0])*0.01
                for dtw,ymin,ymax in zip(track_dtws, ys[:-1], ys[1:]):
                    self._plot_bases(fig, dtw, ymin+dy, ymax-dy, row, col)

            colors = self.prms.track_colors
            dtw_kws = [{"legendgroup" : t.name, "showlegend" : False} for t in self.tracks.alns]

        fig.add_trace(go.Scattergl(
            x=samples, y=signal,
            hoverinfo="skip",
            name="Raw Signal",
            mode="markers",
            marker={"size":2, "color":"black"},
            legendrank=0
        ), row=row, col=col)


        if not self.prms.no_model:
            for dtw,color,kw in zip(track_dtws, colors, dtw_kws):
                dtw = dtw.sort_values("start")
                fig.add_trace(go.Scattergl(
                    name = "Model Current",
                    mode = "lines",
                    x=dtw["start"], y=dtw["model_current"],
                    line={"color":color, "width":2, "shape" : "hv"},
                    **kw
                ), row=row, col=col)

        fig.update_yaxes(title_text="Current (pA)", row=row, col=col)
        fig.update_yaxes(fixedrange=self.prms.yaxis_fixed, row=row, col=col)

        return fig

OPTS = (
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt("input", "tracks", nargs="+"),
    Opt(("-o", "--outfile"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),
    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-n", "--max-reads"), "sigplot"),
)

def main(conf):
    """plot a dotplot"""
    io = Tracks(conf=conf)

    fig = Sigplot(tracks, conf=conf).plot()

    fig.write_html(conf.outfile, config={"scrollZoom" : True, "displayModeBar" : True})
