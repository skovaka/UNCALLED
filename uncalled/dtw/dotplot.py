import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter

from ..sigproc import ProcRead
from ..config import Config, ArgParser, ParamGroup, Opt
from ..fast5 import Fast5Reader, parse_read_ids
from ..pafstats import parse_paf
from .dtw import Track, ref_coords
#from .align import GuidedDTW, BcFast5Aln

class DotplotParams(ParamGroup): pass
DotplotParams._def_params(
    ("track", None, str, "DTW aligment track containing reads to plot"),
    ("read_filter", None, list, "List of reads to plot"),
)
Config._EXTRA_GROUPS["dotplot"] = DotplotParams #TODO put in ParamGroup con

OPTS = [
    Opt("track", "dotplot"),
    Opt(("-R", "--ref-bounds"), "align", type=ref_coords),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
]

def main(conf):
    """Plot dotplots of alignments from tracks produced by `align` or `convert`"""

    track = Track(conf.dotplot.track, 'r')
    conf = track.conf

    fast5s = Fast5Reader(conf=conf)

    ref_bounds = conf.align.ref_bounds

    read_filter = set(conf.fast5_reader.read_filter)

    mm2s = {p.qr_name : p
             for p in parse_paf(
                conf.align.mm2_paf,
                ref_bounds,
                read_filter=read_filter
    )}

    for read_id in track.read_ids:
        if read_id in fast5s:
            fast5_read = fast5s[read_id]
            proc_read = ProcRead(fast5_read, conf=conf)
            bcaln = BcFast5Aln(proc_read, mm2s[read.id])

            aln = track.get_aln(read_id)

            #TODO shouldn't need GuidedDTW, just read_aln
            dtw = GuidedDTW(
                self.idx, 
                proc_read, 
                bcaln, 
                dtw_events=aln,
                ref_bounds=ref_bounds
            )

class Dotplot:
    def __init__(self, index, dtw, out_prefix=None, cursor=None):
        self._init_plot()

        self.index = index

        self.ax_dot.set_ylabel("%s (%s)" % (dtw.bcaln.ref_name, "+" if dtw.bcaln.is_fwd else "-"), fontsize=12)

        self.ax_dot.yaxis.set_major_formatter(FuncFormatter(self._tick_formatter))
        
        self.ax_sig.set_title(dtw.read.id)

        dtw.bcaln.plot_scatter(self.ax_dot, False, samp_min=dtw.samp_min, samp_max=dtw.samp_max)

        #TODO move to plot_cursor, don't put in constructor
        if cursor is not None:
            cursor_kw = {
                'color' : 'red', 
                'alpha' : 0.5
            }
            if dtw.bcaln.flip_ref:
                cursor_ref = np.abs(dtw.bcaln.y_min + cursor)
            else:
                cursor_ref = cursor - dtw.bcaln.y_min
            i = dtw.dtw['ref'].searchsorted(cursor_ref)
            cursor_samp = dtw.dtw.iloc[i]['sample'] + dtw.dtw.iloc[i]['length']/2
            
            self.ax_dot.axvline(cursor_samp, **cursor_kw),
            self.ax_dot.axhline(cursor_ref,  **cursor_kw)

        dtw.plot_dotplot(self.ax_dot)

        dtw.plot_dtw_events(self.ax_sig, self.ax_padiff)

        self.fig.tight_layout()

    def _tick_formatter(self, x, pos):
        return self.index.mirror_to_ref(int(x))

    def _init_plot(self):
        matplotlib.use("TkAgg")
        plt.style.use(['seaborn'])

        self.fig = plt.figure()
        gspec = self.fig.add_gridspec(
                ncols=2, nrows=2, 
                width_ratios=[3,1],
                height_ratios=[1,3] 
        )

        self.ax_sig = self.fig.add_subplot(gspec[0,0])
        self.ax_dot = self.fig.add_subplot(gspec[1,0])
        self.ax_padiff = self.fig.add_subplot(gspec[1,1])

        def sharex(ax1, ax2):
            ax1.get_shared_x_axes().join(ax1,ax2)

        def sharey(ax1, ax2):
            ax1.get_shared_y_axes().join(ax1,ax2)

        sharex(self.ax_dot, self.ax_sig)
        sharey(self.ax_dot, self.ax_padiff)

        fontsize=12
        self.ax_sig.set_ylabel("Current (pA)", fontsize=fontsize)
        self.ax_dot.set_xlabel("Raw Sample", fontsize=fontsize)
        self.ax_padiff.set_xlabel("Abs. pA Difference", fontsize=fontsize)

        self.ax_sig.xaxis.set_major_formatter(NullFormatter())
        self.ax_padiff.yaxis.set_major_formatter(NullFormatter())


    def show(self):
        #self._init_plot()

        plt.show()
        plt.close()
    
    def save(self, out_prefix):
        #self._init_plot()

        suffix = read.id + ".png"
        if os.path.isdir(out_prefix):
            out_fname = os.path.join(out_prefix, suffix)
        else:
            out_fname = out_prefix + suffix

        self.fig.set_size_inches(8, 6)
        self.fig.savefig(out_fname, dpi=100)
        plt.close()
