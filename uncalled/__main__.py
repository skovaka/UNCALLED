from time import time 
t = time()

import sys
import numpy as np
import pandas as pd
import textwrap
from .argparse import ArgParser, Opt, MutexOpts, CONFIG_PARAM, FAST5_PARAM, comma_split, ref_coords
from .fast5 import parse_read_ids
from .index import str_to_coord
#from ..argparse import Opt, comma_split, ref_coords

#from . import index, map, realtime, sim, pafstats, dtw
from . import index, pafstats
from .rt import realtime, map, sim
from .dtw import dtw, convert, db
from .vis import browser, dotplot, sigplot, trackplot
from .stats import refstats, layerstats, _readstats
print("LOAD", time() - t)

BWA_OPTS = (
    Opt("bwa_prefix", "mapper"),
    Opt(("-p", "--idx-preset"), "mapper"),
)

FAST5_OPTS = (
    Opt(FAST5_PARAM, "fast5_reader", nargs="+", type=str),
    Opt(("-r", "--recursive"), "fast5_reader", action="store_true"),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-x", "--fast5-index"), "fast5_reader"),
    Opt(("-n", "--max-reads"), "fast5_reader")
)

MAPPER_OPTS = (
    Opt(("-t", "--threads"), ""),
    Opt("--num-channels", "read_buffer"),
    Opt(("-c", "--max-chunks"), "read_buffer"),
    Opt("--chunk-time", "read_buffer"),
    Opt("--config",
        type = str, 
        default = None, 
        required = False, 
        help = "Config file",
        dest = CONFIG_PARAM
    ),
    Opt("--rna", fn="set_r94_rna")
)

RT_OPTS = (
    MutexOpts("realtime_mode", [
        Opt(("-D", "--deplete"), fn="set_rt_deplete"),
        Opt(("-E", "--enrich"), fn="set_rt_enrich"),
    ]),
    MutexOpts("active_chs", [
        Opt("--even", fn="set_active_chs_even", help="world"),
        Opt("--odd", fn="set_active_chs_odd", help="Hello"),
    ])
)

REALTIME_OPTS = BWA_OPTS + MAPPER_OPTS + (
    Opt("--host", "realtime"),
    Opt("--port", "realtime"),
    Opt("--duration", "realtime"),
) + RT_OPTS

CONVERT_OPTS = (
    Opt("index_prefix", "tracks"),
    Opt(FAST5_PARAM, "fast5_reader", nargs="+", type=str),
    Opt(("-r", "--recursive"), "fast5_reader", action="store_true"),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-n", "--max-reads"), "fast5_reader"),

    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("-a", "--append"), "tracks.io", action="store_true"),
    Opt(("-o", "--db-out"), "tracks.io", required=True),
)

NANOPOLISH_OPTS = CONVERT_OPTS + (
    Opt(("-x", "--fast5-index"), "fast5_reader", required=True), 
    Opt("eventalign_tsv", type=str, default=None, help="Nanopolish eventalign output (should include")
)

NEW_OPTS = (
    Opt("index_prefix", "tracks"),
    #MutexOpts("input", [
        Opt("--db-in", "tracks.io"),
        Opt("--eventalign-in", "tracks.io", nargs="?", const="-"),
    #]),

    #MutexOpts("output", [
        Opt("--db-out", "tracks.io"),
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
    #]),

    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("-a", "--append"), "tracks.io", action="store_true"),
)

DB_OPT = Opt("db_in", "tracks.io", help="Track database file")

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
    Opt(("-F", "--fast5-files"), "fast5_reader", type=comma_split),
    Opt(("-r", "--recursive"), "fast5_reader", action="store_true"),
)
MERGE_OPTS = (
    Opt("dbs", nargs="+", type=str, help="Database files to merge. Will write to the first file if \"-o\" is not specified. "),
    #Opt(("-n", "--track_names"), nargs="+", help="Names of tracks to merge. Will merge all tracks if not specified"),
    Opt(("-o", "--db-out"), "tracks.io", type=str, default=None, help="Output database file. Will output to the first input file if not specified"),
)

COMPARE_OPTS = (
    Opt("db_in", "tracks.io"),
    #Opt(("-l", "--read-filter"), "tracks", nargs="+", type=str),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-R", "--ref-bounds"), "tracks"),
    Opt(("-b", "--bcaln"), action="store_true", help="Compare against basecalled alignment. If two tracks input will look for \"bcaln\" group in second track, otherwise will look in the first track."),
    Opt(("-s", "--save"), action="store_true", help="Will save in database if included, otherwise will output to TSV. Will be associated with the first track listed."),
    Opt(("-j", "--jaccard"), action="store_true", help="Will compute per-reference raw sample jaccard distances. Output by default if no other statistics are specified."),
    Opt(("-d", "--mean-ref-dist"), action="store_true", help="Will compute mean reference coordinate distances between raw samples of alignments of the same read. Output by default if no other statistics are specified."),
    #Opt(("-o", "--output"), choices=["db", "tsv"], help="If \"db\" will output into the track database. If \"tsv\" will output a tab-delimited file to stdout."),
)

DUMP_OPTS = (
    Opt("db_in", "tracks.io", type=str),
    Opt("layers", nargs="+",  help="Layers to retrieve or compute"),
    Opt(("-R", "--ref-bounds"), "tracks", type=ref_coords),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
)

CMDS = {
    "index" : ("Build an index from a FASTA reference", (
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
     )), 

    "map" : ("Rapidly map fast5 read signal to a reference",
        BWA_OPTS + FAST5_OPTS + MAPPER_OPTS
    ), 

    "sim" : ("Simulate real-time targeted sequencing",
        BWA_OPTS + (
        Opt("fast5s", 
            nargs = '+', 
            type = str, 
            help = "Reads to unc. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
        ),
        Opt(("-r", "--recursive"), 
            action = "store_true"
        ),
        Opt("--ctl-seqsum", "simulator"),
        Opt("--unc-seqsum", "simulator"),
        Opt("--unc-paf",    "simulator"),
        Opt("--sim-speed",  "simulator"),
        ) + MAPPER_OPTS
    ), 

    "pafstats" : ("Estimate speed and accuracy from an Uncalled PAF file", (
        Opt("infile",  
            type = str, 
            help = "PAF file output by UNCALLED"
        ),
        Opt(("-n", "--max-reads"), 
            required = False, 
            type = int, 
            default = None, 
            help = "Will only look at first n reads if specified"
        ),
        Opt(("-r", "--ref-paf"), 
            required = False, 
            type = str, 
            default = None, 
            help = "Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA."
        ),
        Opt(("-a", "--annotate"), 
            action = 'store_true', 
            help = "Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag"
        ),
    )), 
    "dtw" : ("Perform DTW alignment guided by basecalled alignments", (
        Opt("index_prefix", "tracks"),) + FAST5_OPTS + (
        Opt(("-m", "--mm2-paf"), "dtw", required=True),
        Opt(("-o", "--db-out"), "tracks.io"),
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
        Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
        Opt(("-a", "--append"), "tracks.io", action="store_true"),
        Opt("--bc-cmp", action="store_true", help="Compute distance from basecalled alignment and store in database"),
        Opt(("-p", "--pore-model"), "pore_model", "name"),
        Opt("--full-overlap", "tracks", action="store_true"),
        #Opt(("-S", "--mask-skips"), "dtw", action="store_true"),
        Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
        Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
        #Opt("--method", "dtw", choices=METHODS.keys()),
        Opt(("-i", "--iterations"), "dtw"),
        Opt(("-c", "--cost-fn"), "dtw", choices=["abs_diff","z_score","norm_pdf"]),
        Opt("--skip-cost", "dtw"),
        Opt("--stay-cost", "dtw"),
        Opt("--move-cost", "dtw"),
        Opt(("-b", "--band-width"), "dtw"),
        Opt(("-s", "--band-shift"), "dtw"),
        Opt(("-N", "--norm-len"), "normalizer", "len", default=0),
    )), 
    "convert" : ("Import DTW alignments produced by other tools into a database", {
        "new" : ("New convert interface", NEW_OPTS), 
        "nanopolish" : ("Convert from nanopolish eventalign TSV to uncalled DTW track", NANOPOLISH_OPTS), 
        "tombo" : ("Convert from Tombo resquiggled fast5s to uncalled DTW track", CONVERT_OPTS)
    }), 
    "db" : ("Edit, merge, and ls alignment databases", {
        "ls"     : ("", LS_OPTS), 
        "delete" : ("", DELETE_OPTS), 
        "merge"  : ("", MERGE_OPTS), 
        "edit"   : ("", EDIT_OPTS), 
    }),
    "refstats" : ("Calculate per-reference-coordinate statistics", (
        Opt("db_in", "tracks.io", type=str),
        Opt("layers", "tracks", type=comma_split,
            help="Comma-separated list of layers over which to compute summary statistics"),# {%s}" % ",".join(LAYERS.keys())),
        Opt("refstats", type=comma_split,
            help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_REFSTATS)),
        Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
        Opt(("-C", "--min-coverage"), "tracks"),
        Opt(("--ref-chunksize"), "tracks.io"),
        Opt(("-c", "--cov"), action="store_true", help="Output track coverage for each reference position"),
        Opt(("-v", "--verbose-refs"), action="store_true", help="Output reference name and strand"),
    )),
    "layerstats" : ("Compute, compare, and query alignment layers", {
        "compare" : ("Compute distance between alignments of the same reads", COMPARE_OPTS),
        "dump" : ("Output DTW alignment paths and statistics", DUMP_OPTS),
    }),
    "dotplot" : ("Plot signal-to-reference alignment dotplots", (
        Opt("db_in", "tracks.io"),
        Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot.", required=True),
        Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
        Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
        Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
        Opt(("-L", "--layers"), "dotplot", "layers", type=comma_split),
        Opt(("-b", "--bcaln-track"), "dotplot"),
        Opt(("-p", "--pore-model"), "pore_model", "name"),
        Opt(("--multi-background"), "sigplot", action="store_true"),
        Opt(("--show-events"), "sigplot", action="store_true"),
        Opt(("--no-model"), "sigplot", action="store_true"),
        Opt(("--bcaln-error", "-e"), "dotplot", action="store_true"),
    )),
    "trackplot" : ("Plot alignment tracks and per-reference statistics", (
        Opt("db_in", "tracks.io"),
        Opt("ref_bounds", "tracks", type=str_to_coord),
        Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
        Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
        Opt(("-H", "--panel-heights"), "trackplot", nargs="+", type=int),
        Opt(("--shared-refs-only"), "tracks", action="store_true"),
        Opt(("--shared-reads-only"), "tracks", action="store_true"),
        Opt(("--share-reads"), "trackplot", action="store_true"),
        Opt(("--hover-read"), "trackplot", action="store_true"),

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
        Opt(("-o", "--outfile"), "trackplot"),
    )),
    "browser" : ("Interactive signal alignment genome browser", (
        Opt("db_in", "tracks.io"),
        Opt("ref_bounds", "tracks", type=str_to_coord),
        #Opt("layer", "trackplot", default="current", nargs="?"),
        Opt(("-r", "--refstats"), "tracks", default=None, type=comma_split),
        Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
        Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
        Opt("--pore-model", "pore_model", "name"),
        Opt(("-p", "--browser-port"), help="Browser port", default=8000),
        Opt(("-o", "--outfile"), "trackplot"),
    )),
}

SUBCMDS = [
    index, 
    map, sim, pafstats, #realtime, 
    dtw, convert, db,
    browser, dotplot, sigplot, trackplot,
    refstats, _readstats, layerstats,
]

_help_lines = [
    "Utility for Nanopore Current ALignment to Large Expanses of DNA", "",
    "subcommand options:",
    "General:",
    "\tindex      Build an index from a FASTA reference",
    "Real-Time Enrichment (Rapid Signal Mapping):",
#    "\trealtime   " + realtime.main.__doc__,
    "\tmap        Rapidly map fast5 read signal to a reference",
    "\tsim        Simulate real-time targeted sequencing",
    "\tpafstats   Estimate speed and accuracy from an Uncalled PAF file",
    "Dynamic Time Warping (DTW) Alignment:",
    "\tdtw        Perform DTW alignment guided by basecalled alignments",
    "\tconvert    Import DTW alignments produced by other tools into a database",
    "\tdb         " + db.__doc__.split("\n")[0], "",
    "DTW Analysis:",
    "\trefstats   " + refstats.main.__doc__,
    #"\treadstats  " + _readstats.main.__doc__,
    "\tlayerstats " + layerstats.__doc__.split("\n")[0],"",
    "DTW Visualization:",
    "\tdotplot    " + dotplot.main.__doc__,
    "\ttrackplot  " + trackplot.main.__doc__,
    "\tbrowser    " + browser.main.__doc__,
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
    t = time()
    parser = ArgParser(SUBCMDS, HELP)
    print("Parser", time()-t)
    t = time()

    cmd, conf = parser.parse_args()

    print("Parsed", time()-t)
    t = time()

    if cmd is not None:
        ret = cmd(conf=conf)

        if isinstance(ret, pd.DataFrame):
            ret.to_csv(sys.stdout, sep="\t")

    else:
        parser.print_help()

    print("done", time()-t)
    t = time()
if __name__ == "__main__":
    main()
