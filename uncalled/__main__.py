import sys
import textwrap
import importlib
import os

from .argparse import ArgParser, Opt, MutexOpts, CONFIG_PARAM, FAST5_PARAM, comma_split, ref_coords

from .index import RefCoord, str_to_coord
import pandas as pd

import cProfile

def parse_read_ids(reads):
    if reads is None:
        return []

    if isinstance(reads, str):
        if os.path.exists(reads):
            with open(reads) as reads_in:
                return [line.split()[0] for line in reads_in]
        else:
            return reads.split(",")

    return list(reads)

CONFIG_OPT = Opt(("-C", "--config"), type=str, default=None, required=False, help="Configuration file in TOML format", dest = CONFIG_PARAM)

FAST5_OPTS = (
    Opt(FAST5_PARAM, "read_index", nargs="+", type=str),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-n", "--max-reads"), "read_index")
)

DTW_OPTS = (
    Opt(("-p", "--processes"), "tracks.io"),
    Opt("--bam-chunksize", "tracks.io"),
    Opt("--bam-in", "tracks.io", nargs="?", const="-", required=True),
    Opt(("--out-name", "-o"), "tracks.io"),
    Opt("ref_index", "tracks"), #+ FAST5_OPTS + (
    
    Opt("--flowcell", "pore_model"),
    Opt("--kit", "pore_model"),

    Opt(FAST5_PARAM, "read_index", nargs="+", type=str),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt(("-l", "--read-filter"), "read_index"),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-n", "--max-reads"), "tracks"),

    Opt("--count-events", "tracks", action="store_true"),
    Opt("--del-max", "dtw"),
    Opt("--ins-max", "dtw"),
    Opt("--mask-skips", "tracks", nargs="?", const="all"),
    Opt("--mask-indels", "tracks"),

    #Opt("eventalign_tsv", type=str, default=None, help="Nanopolish eventalign output (should include"),
    Opt("--ordered-out", "tracks.io", action="store_true"),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("--kmer-shift"), "pore_model", "shift", default=None),
    Opt("--save-bands", "dtw", action="store_true"),
    Opt("--full-overlap", "tracks", action="store_true"),
    #Opt(("-S", "--mask-skips"), "dtw", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    #Opt("--method", "dtw", choices=METHODS.keys()),
    Opt(("--iterations"), "dtw"),
    Opt(("-c", "--cost-fn"), "dtw", choices=["abs_diff","z_score","norm_pdf"]),
    Opt("--skip-cost", "dtw"),
    Opt("--stay-cost", "dtw"),
    Opt("--move-cost", "dtw"),
    Opt(("-b", "--band-width"), "dtw"),
    Opt(("-s", "--band-shift"), "dtw"),
    Opt("--mvcmp-mask", "tracks"),
    Opt("--max-norm-dist", "tracks"),
    Opt("--min-aln-length", "tracks"),
    Opt(("-N", "--norm-mode"), "normalizer", "mode", choices=["ref_mom", "model_mom"]),
    Opt("--norm-median", "normalizer", "median", action="store_true"),
    Opt("--norm-full", "normalizer", "full_read", action="store_true"),
    Opt("--unmask-splice", "dtw", action="store_true"),
    CONFIG_OPT,
)

CONVERT_OPTS = (
    Opt("ref_index", "tracks", nargs="?"),
    Opt(("-p", "--processes"), "tracks.io"),
    MutexOpts("input", [
        Opt("--eventalign-in", "tracks.io", type=comma_split, nargs="?", const="-"),
        Opt("--tombo-in", "tracks.io", type=comma_split, action="extend"),
    ]),
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend"),
    Opt(("-t", "--tracks"), "tracks.io", "input_names", type=comma_split),
    Opt(("-m", "--pore-model"), "pore_model", "name"),
    Opt("--kmer-shift", "pore_model", "shift"),
    Opt("--bam-chunksize", "tracks.io"),

    MutexOpts("output", [
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
        Opt("--tsv-out", "tracks.io", nargs="?", const="-"),
        Opt("--bam-out", "tracks.io", nargs="?", const="-"),
    ]),
    Opt(("--out-name", "-o"), "tracks.io"),

    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtw"),
    Opt("--tsv-na", "tracks.io", nargs="?", const="-"),
    Opt("--tsv-noref", "tracks.io", action="store_true"),

    Opt("--eventalign-flags", "tracks.io", type=comma_split),
    Opt("--mask-skips", "tracks", nargs="?", const="all"),

    Opt("--flowcell", "pore_model"),
    Opt("--kit", "pore_model"),
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt(("-l", "--read-filter"), "tracks"),
    Opt(("-x", "--read-index"), "read_index", required=False),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("-a", "--append"), "tracks.io", action="store_true"),
    CONFIG_OPT,
)

COMPARE_OPTS = (
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend"),
    Opt(("-t", "--tracks"), "tracks.io", "input_names", type=comma_split),
    #Opt(("-l", "--read-filter"), "tracks", nargs="+", type=str),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-R", "--ref-bounds"), "tracks"),
    Opt(("-m", "--moves"), action="store_true", help="Compare against basecalled alignment. If two tracks input will look for \"moves\" group in second track, otherwise will look in the first track."),
    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtwcmp,mvcmp"),
    Opt("--tsv-na", "tracks.io", nargs="?", const="-"),
    Opt("--tsv-noref", "tracks.io", action="store_true"),
    MutexOpts("output", [
        Opt("--tsv-out", "tracks.io", nargs="?", const="-"),
        Opt("--bam-out", "tracks.io", nargs="?", const="-"),
    ]),
    #Opt(("-s", "--save"), action="store_true", help="Will save in database if included, otherwise will output to TSV. Will be associated with the first track listed."),
    #Opt(("-j", "--jaccard"), action="store_true", help="Will compute per-reference raw sample jaccard distances. Output by default if no other statistics are specified."),
    #Opt(("-d", "--mean-ref-dist"), action="store_true", help="Will compute mean reference coordinate distances between raw samples of alignments of the same read. Output by default if no other statistics are specified."),
    CONFIG_OPT,
    #Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
)

ALL_REFSTATS = {"min", "max", "mean", "median", "stdv", "var", "skew", "kurt", "q25", "q75", "q5", "q95", "KS"}
REFSTATS_OPTS = (
    Opt("layers", "tracks", type=comma_split,
        help="Comma-separated list of layers over which to compute summary statistics"),# {%s}" % ",".join(LAYERS.keys())),
    Opt("refstats", type=comma_split,
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_REFSTATS)),
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend", nargs="?", const="-"), #, required=True
    Opt(("-t", "--tracks"), "tracks.io", "input_names", type=comma_split),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt("--min-coverage", "tracks"),
    Opt(("--ref-chunksize"), "tracks.io"),
    Opt(("--aln-chunksize"), "tracks.io"),
    Opt(("-c", "--cov"), action="store_true", help="Output track coverage for each reference position"),
    Opt("--ref-index", "tracks", "ref_index"), 
    Opt(("-m", "--pore-model"), "pore_model", "name", default=None),
    Opt(("-p", "--processes"), "tracks.io"),
)

DTW_CMD_OPTS = DTW_OPTS + (
    MutexOpts("output", [
        Opt("--tsv-out", "tracks.io", nargs="?", const="-"),
        Opt("--bam-out", "tracks.io", nargs="?", const="-"),
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
    ]),
    Opt(("-m", "--pore-model"), "pore_model", "name", default=None),
    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtw"),
    Opt("--tsv-na", "tracks.io", nargs="?", const="-"),
    Opt("--tsv-noref", "tracks.io", action="store_true"),
    Opt("--eventalign-flags", "tracks.io", type=comma_split),
    Opt("--mvcmp", action="store_true", help="Compute distance from basecalled alignment and store in database"),
)

TRAIN_OPTS = (
    Opt(("-i", "--train-iterations"), "train", "iterations"), 
    Opt(("-m", "--init-model"), "train"),
    Opt("--init-mode", "train"),
    Opt(("-k", "--kmer-len"), "train"),
    Opt("--kmer-samples", "train"), 
    Opt("--buffer-size", "train"), 
    Opt(("-d", "--max-moves-dist"), "train"), 
    Opt(("--use-median"), "train", action="store_true"), 
    Opt("--out-dir", "tracks.io", "model_dir"),
    Opt(("-a", "--append"), "train", action="store_true"),
    Opt("--skip-dtw", "train", action="store_true"),
) + DTW_OPTS


READSTATS_OPTS = (
    Opt("stats", "readstats", type=comma_split),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    #Opt(("-p", "--pca-components"), "readstats"),
    #Opt(("-L", "--pca-layer"), "readstats"),
    Opt(("-s", "--summary-stats"), "readstats", type=comma_split),
    CONFIG_OPT,
)

REFPLOT_OPTS = (
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend"),
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-L", "--layer"), "refplot"),
    CONFIG_OPT,
    #Opt(("-o", "--outfile"), "vis"),
)

DOTPLOT_OPTS = (
    MutexOpts("input", [
        Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
        Opt("--bam-in", "tracks.io", nargs="?", const="-", type=comma_split, action="extend"),
        Opt("--eventalign-in", "tracks.io", nargs="?", const="-", type=comma_split, action="extend"),
    ]),

    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),

    Opt("--ref", "tracks", "ref_index"), 
    Opt("--names", "tracks.io", "input_names", type=comma_split), 
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),

    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    #Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", "layers", type=comma_split),
    Opt(("-b", "--moves-track"), "dotplot"),
    Opt(("-p", "--pore-model"), "pore_model", "name", default=None),
    Opt(("--multi-background"), "sigplot", action="store_true"),
    Opt(("--show-events"), "sigplot", action="store_true"),
    Opt(("--show-bands"), "dotplot", action="store_true"),
    Opt(("--no-model"), "sigplot", action="store_true"),
    Opt(("--moves-error", "-e"), "dotplot", action="store_true"),
    CONFIG_OPT,
)

TRACKPLOT_OPTS = (
    Opt("ref_bounds", "tracks", type=str_to_coord),
    MutexOpts("input", [
		Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
		Opt("--bam-in", "tracks.io", nargs="?", const="-", type=comma_split, action="extend"),
	]),

    Opt("--ref", "tracks", "ref_index"), 
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt("--pore-model", "pore_model", "name"),

    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    #Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-H", "--panel-heights"), "trackplot", nargs="+", type=int),
    Opt(("--shared-refs-only"), "tracks", action="store_true"),
    Opt(("--shared-reads-only"), "tracks", action="store_true"),
    Opt(("--share-reads"), "trackplot", action="store_true"),
    Opt(("--hover-read"), "trackplot", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
    CONFIG_OPT,
)

BROWSER_OPTS = (
    Opt("ref_bounds", "tracks", type=RefCoord),
    #MutexOpts("input", [
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend", nargs="?", const="-"),
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend", nargs="?", const="-"), #, required=True
	#]),

    Opt("--ref-index", "tracks", "ref_index"), 
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),

    #Opt("layer", "trackplot", default="current", nargs="?"),
    #Opt(("-r", "--refstats"), "tracks", default=None, type=comma_split),
    Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt("--pore-model", "pore_model", "name"),
    Opt(("-p", "--port"), help="Browser port", default=8000),
    Opt(("-o", "--outfile"), "trackplot"),
    CONFIG_OPT,
)

def panel_opt(name):
    return (lambda arg: (name, arg))

TRACKPLOT_PANEL_OPTS = (
    Opt("--mat", dest="panels",
        metavar="LAYER", action="append", type=panel_opt("mat"),
        help="Display a ref-by-read matrix of specified alignment layer"), 

    Opt("--box", dest="panels", #"trackplot", "panels", 
        metavar="LAYER", action="append", type=panel_opt("box"),
        help="Display a boxplot of specified layer"), 

    Opt("--line", dest="panels", #"trackplot", "panels", 
        metavar="LAYER.STAT", action="append", type=panel_opt("line"),
        help="Display a line plot of specifed layer summary statistic"), 

    Opt("--scatter", dest="panels", #"trackplot", "panels", 
        metavar="LAYER.STAT", action="append", type=panel_opt("scatter"),
        help="Display a line plot of specifed layer summary statistic"), 
)


CMDS = {
    "dtw" : ("dtw.dtw", 
        "Perform DTW alignment guided by basecalled alignments", DTW_CMD_OPTS), 
    "convert" : ("dtw.io", "Convert between signal alignment file formats", CONVERT_OPTS),
    "train" : ("dtw.train", 
        "Iteratively train a new k-mer pore model", TRAIN_OPTS), 
    "refstats" : ("stats.refstats", "Calculate per-reference-coordinate statistics", REFSTATS_OPTS),
    "readstats" : ("stats.readstats", "", READSTATS_OPTS),
    "compare" : ("stats.layerstats", "Compute distance between alignments of the same reads", COMPARE_OPTS),
    "dotplot" : ("vis.dotplot", "Plot signal-to-reference alignment dotplots", DOTPLOT_OPTS),
    "refplot" : ("vis.refplot", "Plot alignment tracks and per-reference statistics", REFPLOT_OPTS),
    "trackplot" : ("vis.trackplot", "Plot alignment tracks and per-reference statistics", TRACKPLOT_OPTS+TRACKPLOT_PANEL_OPTS),
    "browser" : ("vis.browser", "Interactive signal alignment genome browser", BROWSER_OPTS),
}

_help_lines = [
    "Utility for Nanopore Current ALignment to Large Expanses of DNA", "",
    "subcommand options:",
    "Dynamic Time Warping (DTW) Alignment:",
    "\tdtw        Perform DTW alignment guided by basecalled alignments",
    "\ttrain      Train new k-mer pore models",
    "\tconvert    Convert between signal alignment file formats",
    "DTW Analysis:",
    "\trefstats   Calculate per-reference-coordinate statistics",
    "\treadstats  Perform per-read analyses of DTW alignments",
    #"\tlayerstats Compute, compare, and query alignment layers", "",
    "DTW Visualization:",
    "\tdotplot    Plot signal-to-reference alignment dotplots",
    "\ttrackplot  Plot alignment tracks and per-reference statistics",
    "\tbrowser    Interactive signal alignment genome browser",
    #"\tsigplot    " + sigplot.main.__doc__,
]

HELP = "\n".join([
    textwrap.fill(line,
        width=75,
        drop_whitespace=False, 
        replace_whitespace=False, 
        tabsize=2, 
        subsequent_indent=" "*13)
    for line in _help_lines])


def main():
    #parser = ArgParser(SUBCMDS, HELP)
    parser = ArgParser(CMDS, HELP)

    module, cmd, conf = parser.parse_args()


    if module is not None:

        m = importlib.import_module(f".{module}", "uncalled")
        fn = getattr(m, cmd)
        if conf.cprof is not None:
            cProfile.runctx(f"fn(conf=conf)",
             {"fn" : fn, "conf": conf},{},
             conf.cprof)
        else:
            ret = fn(conf=conf)
            if isinstance(ret, pd.DataFrame):
                ret.to_csv(sys.stdout, sep="\t")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
