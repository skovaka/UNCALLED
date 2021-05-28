import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter

def dotplot(dtw, out_prefix=None, cursor=None):

    matplotlib.use("TkAgg")

    plt.style.use(['seaborn'])

    t = time.time()

    #dims = (4,8)
    #fig = plt.figure()
    #ax_sig =  plt.subplot2grid(dims, (0,0), colspan=6, rowspan=1)
    #ax_dot =  plt.subplot2grid(dims, (1,0), colspan=6, rowspan=3)
    #ax_padiff = plt.subplot2grid(dims, (1,6), colspan=1, rowspan=3)
    #ax_bcdist = plt.subplot2grid(dims, (1,7), colspan=1, rowspan=3)

    dims = (4,4)
    fig = plt.figure()
    ax_sig =  plt.subplot2grid(dims, (0,0), colspan=3, rowspan=1)
    ax_dot =  plt.subplot2grid(dims, (1,0), colspan=3, rowspan=3)
    ax_padiff = plt.subplot2grid(dims, (1,3), colspan=1, rowspan=3)


    def sharex(ax1, ax2):
        ax1.get_shared_x_axes().join(ax1,ax2)

    def sharey(ax1, ax2):
        ax1.get_shared_y_axes().join(ax1,ax2)

    sharex(ax_dot, ax_sig)
    sharey(ax_dot, ax_padiff)


    fontsize=12
    ax_sig.set_ylabel("Current (pA)", fontsize=fontsize)
    ax_dot.set_xlabel("Raw Sample", fontsize=fontsize)
    ax_dot.set_ylabel("%s (%s)" % (dtw.bcaln.rf_name, "+" if dtw.bcaln.is_fwd else "-"), fontsize=fontsize)
    ax_padiff.set_xlabel("Abs. pA Difference", fontsize=fontsize)

    #ax_bcdist.set_xlabel("BC Distance", fontsize=fontsize)

    #ax_sig.set_xticks([])
    ax_sig.xaxis.set_major_formatter(NullFormatter())
    ax_padiff.yaxis.set_major_formatter(NullFormatter())
    #ax_bcdist.yaxis.set_major_formatter(NullFormatter())

    ax_dot.yaxis.set_major_formatter(FuncFormatter(dtw.bcaln.ref_tick_fmt))
    
    ax_sig.set_title(dtw.read.id)
    #dtw.bcaln.plot_step(ax_dot, False, samp_min=dtw.samp_min, samp_max=dtw.samp_max)
    dtw.bcaln.plot_scatter(ax_dot, False, samp_min=dtw.samp_min, samp_max=dtw.samp_max)

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
        
        ax_dot.axvline(cursor_samp, **cursor_kw),
        ax_dot.axhline(cursor_ref,  **cursor_kw)

    dtw.plot_dotplot(ax_dot)

    dtw.plot_dtw_events(ax_sig, ax_padiff)


    if out_prefix is None:
        #mng = plt.get_current_fig_manager()
        #mng.window.maximize()
        plt.show()
        plt.close()

    else:
        suffix = read.id + ".png"
        if os.path.isdir(out_prefix):
            out_fname = os.path.join(out_prefix, suffix)
        else:
            out_fname = out_prefix + suffix

        #fig = plt.gcf()

        fig.set_size_inches(8, 6)
        fig.tight_layout()
        fig.savefig(out_fname, dpi=100)
        plt.close()
