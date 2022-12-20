import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import os

from .. import config
from ..index import str_to_coord
from ..dtw.tracks import Tracks
from ..argparse import Opt, comma_split
from ..fast5 import parse_read_ids
from ..signal_processor import SignalProcessor

class Sigplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("sigplot", *args, **kwargs)
        self.tracks = self.prms.tracks #TODO do this better

        self.sigproc = SignalProcessor(self.tracks.alns[0].model, self.conf)

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

            
    def _plot_bases(self, fig, dtw, model, ymin, ymax, row, col):
        bases = model.kmer_base(dtw["kmer"], 2)
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
                name=model.base_to_char(base),
                legendrank=1
            ), row=row, col=col)
            self._legend.add(color)

    def plot_read(self, fig, read_id, row=1, col=1):

        track_dtws  = list()
        
        current_min = samp_min = np.inf
        current_max = samp_max = 0
        fast5_file = None
        for i,track in enumerate(self.tracks.alns):
            #track_color = self.prms.track_colors[i]
            colors = self.conf.vis.track_colors[i]

            alns = track.alignments.query("@read_id == read_id")
            aln_ids = alns.index


            dtws = list()
            
            for aln_id,aln in alns.iterrows():
                if "fast5" in aln.index:
                    fast5_file = aln.loc["fast5"]

                dtw = track.layers.loc[(slice(None),aln_id),"dtw"].droplevel("aln_id")
                dtw["model_current"] = track.model[dtw["kmer"]]
                dtws.append(dtw)

                max_i = dtw["start"].argmax()
                samp_min = min(samp_min, dtw["start"].min())
                samp_max = max(samp_max, dtw["start"].iloc[max_i] + dtw["length"].iloc[max_i])

                current_min = min(current_min, dtw["model_current"].min())
                current_max = max(current_max, dtw["model_current"].max())

            if len(dtws) > 0:
                track_dtws.append(pd.concat(dtws).sort_index())

        if fast5_file is None:
            fast5 = self.tracks.fast5s[read_id]
        else:
            fast5 = self.tracks.fast5s.get_read(read_id, os.path.basename(fast5_file))

        read = self.sigproc.process(fast5, True)
        signal = read.get_norm_signal(samp_min, samp_max)

        sig_med = np.median(signal)
        sig_win = signal.std()*3 #, signal.min(), signal.max())

        #TODO set global signal min/max
        mask = ((signal >= sig_med - sig_win) &
                (signal <= sig_med + sig_win))
        
        samples = np.arange(samp_min, samp_max)[mask]
        signal = signal[mask]

        current_min = min(current_min, signal.min())
        current_max = max(current_max, signal.max())

        if len(self.tracks) == 1:
            self._plot_bases(fig, track_dtws[0], track.model, current_min, current_max, row, col)
            colors = ["white"]
            dtw_kws = [{}]
        else:

            if self.prms.multi_background:
                ys = np.linspace(current_min, current_max, len(self.tracks)+1)
                dy = (ys[1]-ys[0])*0.01
                for dtw,ymin,ymax in zip(track_dtws, ys[:-1], ys[1:]):
                    self._plot_bases(fig, dtw, track.model, ymin+dy, ymax-dy, row, col)

            colors = self.conf.vis.track_colors
            dtw_kws = [{"legendgroup" : t.name, "showlegend" : False} for t in self.tracks.alns]

        fig.add_trace(go.Scattergl(
            x=samples, y=signal,
            hoverinfo="skip",
            name="Raw Signal",
            mode="markers",
            marker={"size":3, "color":"black", "opacity" : 0.75},
            legendrank=0
        ), row=row, col=col)

        if self.prms.show_events:
            fig.add_trace(go.Scattergl(
                name = "Event Means",
                mode = "lines",
                x=read.df["start"], y=read.df["norm_sig"],
                line={"color":"black", "width":2, "shape" : "hv"},
            ), row=row, col=col)

        if not self.prms.no_model:
            for dtw,color,kw in zip(track_dtws, colors, dtw_kws):
                dtw = dtw.sort_values("start")
                fig.add_trace(go.Scattergl(
                    name = "Model Current",
                    mode = "lines",
                    x=dtw["start"], y=dtw["model_current"],
                    line={"color":color, "width":2.5, "shape" : "hv"},
                    **kw
                ), row=row, col=col)

        fig.update_yaxes(title_text="Current (pA)", row=row, col=col)
        fig.update_yaxes(fixedrange=self.prms.yaxis_fixed, row=row, col=col)

        return fig

OPTS = (
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt("db_in", "tracks.io"),
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
