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
from ..signal_processor import SignalProcessor

class Sigplot:
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("sigplot", *args, **kwargs)
        self.tracks = self.prms.tracks #TODO do this better

        self.conf.normalizer.tgt_mean = 0
        self.conf.normalizer.tgt_stdv = 1
        self.sigproc = SignalProcessor(self.tracks.model, self.conf)

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
        bases = model.kmer_base(dtw["seq","kmer"], 2)
        for base, color in enumerate(self.conf.vis.base_colors):
            base_dtw = dtw[bases == base]
            starts = base_dtw["dtw", "start_sec"]
            ends = starts + base_dtw["dtw", "length_sec"] - (1.0 / model.sample_rate)
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

    #def plot_signal(self, fig, read_id, row, col):

    def plot_read(self, fig, read_id, row=1, col=1):

        track_dtws  = list()
        
        current_min = samp_min = np.inf
        current_max = samp_max = 0
        active_tracks = self.tracks.alignments.index.unique("track")
        for i,track in enumerate(self.tracks.alns):
            #track_color = self.prms.track_colors[i]
            colors = self.conf.vis.track_colors[i]

            if track.name not in active_tracks: continue

            alns = self.tracks.alignments.loc[track.name].query("@read_id == read_id")
            aln_ids = alns.index

            dtws = list()
            seqs = list()
            
            for aln_id,aln in alns.iterrows():
                #dtw = track.layers.loc[(slice(None),aln_id),"dtw"].droplevel("aln.id")
                #seq = track.layers.loc[(slice(None),aln_id),"seq"].droplevel("aln.id")
                layers = self.tracks.layers.loc[track.name].loc[aln_id].sort_values(("dtw","start"))
                dtw = layers["dtw"].dropna()#.droplevel("aln.id")
                seq = layers["seq"].dropna()

                #dtws.append(dtw)
                seqs.append(layers)

                samp_min = min(samp_min, dtw["start"].min())
                samp_max = max(samp_max, dtw["start"].iloc[-1] + dtw["length"].iloc[-1])

                current_min = min(current_min, (dtw["current"]-dtw["current_sd"]*2).min())
                current_max = max(current_max, (dtw["current"]+dtw["current_sd"]*2).max())


            if len(seqs) > 0:
                track_dtws.append(pd.concat(seqs).sort_index())

        read = self.tracks.read_index[read_id]

        if read.empty():        
            pass
            signal = None

        else:

            read = self.sigproc.process(read, True)

            signal = read.get_norm_signal(int(samp_min), int(samp_max))

            sig_med = np.median(signal)
            sig_win = signal.std()*3 #, signal.min(), signal.max())

            #TODO set global signal min/max
            mask = ((signal >= sig_med - sig_win) &
                    (signal <= sig_med + sig_win))
            
            samples = np.arange(samp_min, samp_max)#[mask]
            samp_time = samples / self.tracks.model.sample_rate
            signal = signal#[mask]

            current_min = min(current_min, signal.min())
            current_max = max(current_max, signal.max())

        if len(active_tracks) == 1:
            self._plot_bases(fig, track_dtws[0], self.tracks.model, current_min, current_max, row, col)
            colors = ["white"]
            dtw_kws = [{}]
        else:

            if self.prms.multi_background:
                ys = np.linspace(current_min, current_max, len(self.tracks)+1)
                dy = (ys[1]-ys[0])*0.01
                for dtw,ymin,ymax in zip(track_dtws, ys[:-1], ys[1:]):
                    self._plot_bases(fig, dtw, self.tracks.model, ymin+dy, ymax-dy, row, col)

            colors = self.conf.vis.track_colors
            dtw_kws = [{"legendgroup" : t.name, "showlegend" : False} for t in self.tracks.alns]

        if signal is not None:
            fig.add_trace(go.Scattergl(
                x=samp_time, y=signal,
                hoverinfo="skip",
                name="Raw Signal",
                mode="markers",
                marker={"size":3, "color":"black", "opacity" : 0.75},
                legendrank=0
            ), row=row, col=col)
        else:
            #for dtw,color,kw in zip(track_dtws, colors, dtw_kws):
            for dtw,color in zip(track_dtws, colors):
                dtw = dtw.sort_values(("dtw","start_sec"))
                ymin = dtw["dtw", "current"] - dtw["dtw","current_sd"]
                ymax = dtw["dtw", "current"] + dtw["dtw","current_sd"]
                fig.add_trace(go.Scattergl(
                    name = "Signal (+/-1 stdv)",
                    mode = "lines",
                    fill="tonexty",
                    fillcolor="black",
                    opacity=0.75,
                    legendgroup="signal",
                    x=dtw["dtw", "start_sec"], y=ymin,
                    line={"color":"black", "width":2.5, "shape" : "hv"},
                ), row=row, col=col)

                fig.add_trace(go.Scattergl(
                    name = "Model Current",
                    mode = "lines",
                    x=dtw["dtw", "start_sec"], y=ymax,
                    fill="tonexty",
                    fillcolor="black",
                    legendgroup="signal",
                    showlegend=False,
                    line={"color":"black", "width":2.5, "shape" : "hv"},
                ), row=row, col=col)

        if self.prms.show_events:
            fig.add_trace(go.Scattergl(
                name = "Event Means",
                mode = "lines",
                x=read.df["start_sec"], y=read.df["norm_sig"],
                line={"color":"black", "width":2, "shape" : "hv"},
            ), row=row, col=col)

        if not self.prms.no_model:
            for dtw,color,kw in zip(track_dtws, colors, dtw_kws):
                dtw = dtw.sort_values(("dtw","start_sec"))
                fig.add_trace(go.Scattergl(
                    name = "Model Current",
                    mode = "lines",
                    x=dtw["dtw", "start_sec"], y=dtw["seq", "current"],
                    line={"color":color, "width":2.5, "shape" : "hv"},
                    **kw
                ), row=row, col=col)

        fig.update_yaxes(title_text="Current (norm)", row=row, col=col)
        fig.update_yaxes(fixedrange=self.prms.yaxis_fixed, row=row, col=col)

        return fig

def main(conf):
    """plot a dotplot"""
    io = Tracks(conf=conf)

    fig = Sigplot(tracks, conf=conf).plot()

    fig.write_html(conf.outfile, config={"scrollZoom" : True, "displayModeBar" : True})
