import sys
import textwrap
import importlib

from .argparse import ArgParser, Opt, MutexOpts, CONFIG_PARAM, FAST5_PARAM, comma_split, ref_coords

from .index import str_to_coord
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

INDEX_OPTS = (
        Opt("fasta_filename", "index"),
        Opt(("-o", "--index-prefix"), "index"),
        Opt("--no-bwt", "index", action="store_true"),
        Opt(("-s", "--max-sample-dist"), "index"),
        Opt("--min-samples", "index"),
        Opt("--max-samples", "index"),
        Opt(("-1", "--matchpr1"), "index"),
        Opt(("-2", "--matchpr2"), "index"),
        Opt(("-f", "--pathlen-percentile"), "index"),
        Opt(("-m", "--max-replen"), "index"),
        Opt("--probs", "index"),
        Opt("--speeds", "index"),
        CONFIG_OPT,
     )

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
    Opt("--guppy-in", "tracks.io"),
    Opt("--bam-in", "tracks.io", nargs="?", const="-", required=True),
    Opt(("--out-name", "-o"), "tracks.io"),
    Opt("index_prefix", "tracks"), #+ FAST5_OPTS + (

    Opt(FAST5_PARAM, "read_index", nargs="+", type=str),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt(("-l", "--read-filter"), "tracks"),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-n", "--max-reads"), "tracks"),

    Opt("--del-max", "dtw"),
    Opt("--mask-skips", "tracks", nargs="?", const="all"),
    Opt("--mask-indels", "tracks"),

    #Opt("eventalign_tsv", type=str, default=None, help="Nanopolish eventalign output (should include"),
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
    Opt(("-N", "--norm-mode"), "normalizer", "mode", choices=["ref_mom", "model_mom"]),
    Opt("--norm-median", "normalizer", "median", action="store_true"),
    Opt("--norm-seg", "normalizer", "full_read", action="store_false"),
    CONFIG_OPT,
)

CONVERT_OPTS = (
    Opt("index_prefix", "tracks", nargs="?"),
    MutexOpts("input", [
        Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
        Opt("--bam-in", "tracks.io", type=comma_split, action="extend"),
        Opt("--eventalign-in", "tracks.io", type=comma_split, nargs="?", const="-"),
        Opt("--tombo-in", "tracks.io", type=comma_split, action="extend"),
    ]),
    Opt(("-t", "--tracks"), "tracks.io", "in_names", type=comma_split),

    MutexOpts("output", [
        Opt("--sql-out", "tracks.io"),
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
        Opt("--tsv-out", "tracks.io", nargs="?", const="-"),
        Opt("--bam-out", "tracks.io", nargs="?", const="-"),
    ]),
    Opt(("--out-name", "-o"), "tracks.io"),

    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtw"),

    Opt("--eventalign-flags", "tracks.io", type=comma_split),
    Opt("--mask-skips", "tracks", nargs="?", const="all"),

    Opt("--read-paths", "read_index", "paths", nargs="+", type=str),
    Opt(("-l", "--read-filter"), "tracks"),
    Opt(("-x", "--read-index"), "read_index", required=False),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("-a", "--append"), "tracks.io", action="store_true"),
    CONFIG_OPT,
)

DB_OPT = Opt("sql_in", "tracks.io", help="Track database file")

LS_OPTS = (DB_OPT,)
DELETE_OPTS = (
    DB_OPT,
    Opt("track_name", help="Name of the track to delete"),
)
EDIT_OPTS = (
    DB_OPT,
    Opt("track_name", help="Current track name"),
    Opt(("-N", "--new-name"), default=None, help="New track name"),
    Opt(("-D", "--description"), default=None, help="New track description"),
    Opt(("-F", "--paths"), "read_index", type=comma_split),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
)

MERGE_OPTS = (
    Opt("dbs", nargs="+", type=str, help="Database files to merge. Will write to the first file if \"-o\" is not specified. "),
    #Opt(("-n", "--track_names"), nargs="+", help="Names of tracks to merge. Will merge all tracks if not specified"),
    Opt("--sql-out", "tracks.io", type=str, default=None, help="Output database file. Will output to the first input file if not specified"),
)

COMPARE_OPTS = (
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
    Opt(("-t", "--tracks"), "tracks.io", "in_names", type=comma_split),
    #Opt(("-l", "--read-filter"), "tracks", nargs="+", type=str),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-R", "--ref-bounds"), "tracks"),
    Opt(("-b", "--moves"), action="store_true", help="Compare against basecalled alignment. If two tracks input will look for \"moves\" group in second track, otherwise will look in the first track."),
    Opt(("-s", "--save"), action="store_true", help="Will save in database if included, otherwise will output to TSV. Will be associated with the first track listed."),
    Opt(("-j", "--jaccard"), action="store_true", help="Will compute per-reference raw sample jaccard distances. Output by default if no other statistics are specified."),
    Opt(("-d", "--mean-ref-dist"), action="store_true", help="Will compute mean reference coordinate distances between raw samples of alignments of the same read. Output by default if no other statistics are specified."),
    CONFIG_OPT,
    #Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
)

DUMP_OPTS = (
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
    Opt("layers", nargs="+",  help="Layers to retrieve or compute"),
    Opt(("-R", "--ref-bounds"), "tracks", type=ref_coords),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
)

ALL_REFSTATS = {"min", "max", "mean", "median", "stdv", "var", "skew", "kurt", "q25", "q75", "q5", "q95", "KS"}
REFSTATS_OPTS = (
    Opt("layers", "tracks", type=comma_split,
        help="Comma-separated list of layers over which to compute summary statistics"),# {%s}" % ",".join(LAYERS.keys())),
    Opt("refstats", type=comma_split,
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_REFSTATS)),
    #Opt("--sql-in", "tracks.io"),
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend", nargs="?", const="-"), #, required=True
    Opt(("-t", "--tracks"), "tracks.io", "in_names", type=comma_split),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt("--min-coverage", "tracks"),
    Opt(("--ref-chunksize"), "tracks.io"),
    Opt(("--aln-chunksize"), "tracks.io"),
    Opt(("-c", "--cov"), action="store_true", help="Output track coverage for each reference position"),
)

DTW_CMD_OPTS = DTW_OPTS + (
    MutexOpts("output", [
        Opt("--sql-out", "tracks.io"),
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
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
    Opt("stats", "readstats", type=comma_split),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    #Opt(("-p", "--pca-components"), "readstats"),
    #Opt(("-L", "--pca-layer"), "readstats"),
    Opt(("-s", "--summary-stats"), "readstats", type=comma_split),
    CONFIG_OPT,
)

REFPLOT_OPTS = (
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend"),
    Opt("ref_bounds", "tracks", type=str_to_coord),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
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

    Opt("--ref", "tracks", "index_prefix"), 
    Opt("--read-paths", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),

    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
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

    Opt("--ref", "tracks", "index_prefix"), 
    Opt("--read-paths", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt("--pore-model", "pore_model", "name"),

    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
    Opt(("-H", "--panel-heights"), "trackplot", nargs="+", type=int),
    Opt(("--shared-refs-only"), "tracks", action="store_true"),
    Opt(("--shared-reads-only"), "tracks", action="store_true"),
    Opt(("--share-reads"), "trackplot", action="store_true"),
    Opt(("--hover-read"), "trackplot", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
    CONFIG_OPT,
)

BROWSER_OPTS = (
    Opt("ref_bounds", "tracks", type=str_to_coord),
    #MutexOpts("input", [
    Opt("--sql-in", "tracks.io", type=comma_split, action="extend", nargs="?", const="-"),
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend", nargs="?", const="-"), #, required=True
	#]),

    Opt("--ref", "tracks", "index_prefix"), 
    Opt("--read-paths", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),

    #Opt("layer", "trackplot", default="current", nargs="?"),
    #Opt(("-r", "--refstats"), "tracks", default=None, type=comma_split),
    Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt("--pore-model", "pore_model", "name"),
    Opt(("-p", "--browser-port"), help="Browser port", default=8000),
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
    "index" : ("index", 
        "Build an index from a FASTA reference", INDEX_OPTS), 
    "dtw" : ("dtw.dtw", 
        "Perform DTW alignment guided by basecalled alignments", DTW_CMD_OPTS), 
    "convert" : ("dtw.io", "Convert between signal alignment file formats", CONVERT_OPTS),
    "train" : ("dtw.train", 
        "Iteratively train a new k-mer pore model", TRAIN_OPTS), 
    "db" : (None, 
        """Edit, merge, and ls alignment databases

        subcommand options:
        ls       List all tracks in a database
        delete   Delete a track from a database
        merge    Merge databases into a single file
        edit     Rename, change fast5 paths, or set description""", {

        "ls"     : ("dtw.io.sqlite", "", LS_OPTS), 
        "delete" : ("dtw.io.sqlite", "", DELETE_OPTS), 
        "merge"  : ("dtw.io.sqlite", "", MERGE_OPTS), 
        "edit"   : ("dtw.io.sqlite", "", EDIT_OPTS), 
    }),
    "refstats" : ("stats.refstats", "Calculate per-reference-coordinate statistics", REFSTATS_OPTS),
    "readstats" : ("stats.readstats", "", READSTATS_OPTS),
    "layerstats" : (None, 
            """Compute distance between alignments of the same reads\n"""
            """subcommand options:\n""" 
            """compare  Compute distance between alignments of the same reads\n""", {
        "compare" : ("stats.layerstats", "Compute distance between alignments of the same reads", COMPARE_OPTS),
    }),
    "dotplot" : ("vis.dotplot", "Plot signal-to-reference alignment dotplots", DOTPLOT_OPTS),
    "refplot" : ("vis.refplot", "Plot alignment tracks and per-reference statistics", REFPLOT_OPTS),
    "trackplot" : ("vis.trackplot", "Plot alignment tracks and per-reference statistics", TRACKPLOT_OPTS+TRACKPLOT_PANEL_OPTS),
    "browser" : ("vis.browser", "Interactive signal alignment genome browser", BROWSER_OPTS),
}

_help_lines = [
    "Utility for Nanopore Current ALignment to Large Expanses of DNA", "",
    "subcommand options:",
    "General:",
    "\tindex      Build an index from a FASTA reference",
    "Dynamic Time Warping (DTW) Alignment:",
    "\tdtw        Perform DTW alignment guided by basecalled alignments",
    "\ttrain      Train new k-mer pore models",
    "\tconvert    Convert between signal alignment file formats",
    "\tdb         Edit, merge, and ls alignment databases", "",
    "DTW Analysis:",
    "\trefstats   Calculate per-reference-coordinate statistics",
    "\treadstats  Perform per-read analyses of DTW alignments",
    "\tlayerstats Compute, compare, and query alignment layers", "",
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
