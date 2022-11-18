from . import config
from . import RefCoord

#Index parameter group
class IndexParams(config.ParamGroup):
    _name = "index"
IndexParams._def_params(
    ("fasta_filename", None, str, "FASTA file to index"),
    ("index_prefix", None, str, "Index output prefix. Will use input fasta filename by default"),
    ("no_bwt", False, bool, "Will only generate the pacseq if specified, which is much faster to build. Can only be used with DTW subcommands (NOT map, sim, or realtime)"),
    ("max_sample_dist", 100, int, "Maximum average sampling distance between reference alignments."),
    ("min_samples", 50000, int, "Minimum number of alignments to produce (approximate, due to deterministically random start locations),"),
    ("max_samples", 1000000, int, "Maximum number of alignments to produce (approximate, due to deterministically random start locations),"),
    ("kmer_len", 5, int, "Model k-mer length"),
    ("matchpr1", 0.6334, float, "Minimum event match probability"),
    ("matchpr2", 0.9838, float, "Maximum event match probability"),
    ("pathlen_percentile", 0.05, float, ""),
    ("max_replen", 100, int, ""),
    ("probs", None, str, "Find parameters with specified target probabilites (comma separated)"),
    ("speeds", None, str, "Find parameters with specified speed coefficents (comma separated)"),
)

eventalign_flags = "\", \"".join(["print-read-names", "signal-index", "samples"])#, "scale-events"]

class IOParams(config.ParamGroup):
    _name = "io"
IOParams._def_params(
    #("input", None, None, "Input tracks specifier. Should be in the format <file.db>[:<track1>[,<track2>...]]. If no track names are specified, all tracks will be loaded from the database."),
    #("output", None, None,  "Output track specifier. Should be in the format <file.db>[:<track_name>], where file.db is the output sqlite database. If <track_name> is not specified, will be same as filename (without extension)"),

    #("db_in", None, str, "Input track database"),
    #("db_out", None, str, "Output track database"),

    ("sql_in", None, None, "Input track database"),
    ("sql_out", None, str, "Output track database"),

    ("tsv_out", None, str, "TSV output file (or \"-\"/no argument for stdout)"),
    ("tsv_cols", None, list, "TSV file output alignment layers (comma-separated, can also include \"read_id\""),
    ("tsv_na", "*", str, "Missing value representation for TSV output"),
    ("tsv_noref", False, bool, "Will NOT output reference coordinates to TSV if True"),

    ("bam_in", None, None, "BAM input file (or \"-\"/no argument for stdin)"),
    ("bam_out", None, str, "BAM output file (or \"-\"/no argument for stdout)"),

    ("eventalign_in", None, list, "Eventalign (nanopolish) input file (or \"-\"/no argument for stdin)"),
    ("eventalign_out", None, str, "Eventalign (nanopolish) output file (or \"-\"/no argument for stdout)"),
    ("eventalign_index", None, str, "Nanopolish index file"),
    ("eventalign_flags", [], list, f"Eventalign optional flags (comma-separated list of \"\"{eventalign_flags}\")"),
    ("tombo_in", None, None, "Fast5 files containing Tombo alignments"),

    ("guppy_in", None, str, "Guppy directory containing sequencing summary, BAM files, and FAST5s"),

    ("out_name", None, str, "Output track name (defaults to file basename)"),
    ("in_names", None, list, "Names of tracks to read from input(s)"),

    ("init_track", True, bool, "If true will initialze track"),

    ("output_format", "db", str,  "Output format (db, eventalign)"),
    ("overwrite", False, bool, "Overwrite existing tracks"),
    ("append", False, bool, "Append reads to existing tracks"),
    ("aln_chunksize", 4000, int, "Number of alignments to query for iteration"),
    ("ref_chunksize", 10000, int, "Number of reference coordinates to query for iteration"),
    ignore_toml={"bam_in", "bam_out", "sql_in", "sql_out", "eventalign_in", "eventalign_out", "tombo_in", "eventalign_index", "overwrite", "append"},
    #ignore_toml={"input", "output", "output_format", "overwrite", "append"},
    config_add=False
)

class TracksParams(config.ParamGroup):
    _name = "tracks"
TracksParams._def_params(
    ("io", {}, IOParams, "Track input/output parameters"),
    ("ref_bounds", None, RefCoord, "Only load reads which overlap these coordinates"),
    ("read_filter", None, None, "Only load reads which overlap these coordinates"),
    ("max_reads", None, int, "Only load reads which overlap these coordinates"),

    ("mask_skips", None, None, "Either \"all\" to mask all skips, or \"keep_best\" to mask all but the closest to the model"),
    ("mask_indels", None, int, "Mask positions which overlap basecalled alignment insertions or deletions of this length or longer"),

    ("full_overlap", False, bool, "If true will only include reads which fully cover reference bounds"),
    ("min_coverage", 1, int, "Reference positions with less than this coverage will be excluded from each track (or all tracks if shared_refs_only is true)"),
    ("shared_reads_only", False, bool, "If true will only contain reads shared between all tracks"),
    ("shared_refs_only", False, bool, "If true will only contain reference positions where all tracks have sufficient coverage (see min_coverage)"),

    ("layers", [], None, "Layers to load (e.g. current, dwell, model_diff)"),

    ("load_mat", False, bool, "If true will pivot layers into a matrix"), #TODO change to mat_layers, only do it for them

    ("refstats", None, None, "Per-reference summary statistics to compute for each layer"),
    ("refstats_layers", None, None, "Layers to compute refstats"),

    ("index_prefix", None, str, "BWA index prefix"),
    ("load_fast5s", False, bool, "Load fast5 files"),

    ignore_toml={"ref_bounds", "layers", "full_overlap", "refstats", "refstats_layers", "read_filter"}
)

class VisParams(config.ParamGroup):
    _name = "vis"
VisParams._def_params(
    ("track_colors", ["#AA0DFE", "#1CA71C", "#4676FF", "#d90000"], list, "Track Colors"),
    ("base_colors", ["#80ff80", "#6b93ff", "#ffe543", "#ff8080"], list, "Colors for each base (A,C,G,T/U)"), 
)

class DotplotParams(config.ParamGroup):
    _name = "dotplot"
DotplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("bcaln_track", None, str, "Only display basecalled alignments from this track"),
    ("bcaln_error", False, bool, "Display basecalled alignment errors"),
    ("show_legend", True, bool, "Display legend"),
    ("show_bands", True, bool, "Display DTW bands (if present in DB)"),
    ("select_ref", None, int, "Display a horizontal line at specified reference coordinate"),
    ("layers", [], None, ""),
)

class SigplotParams(config.ParamGroup):
    _name = "sigplot"
SigplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    #("ref_bounds", None, str_to_coord, "DTW aligment tracks"),
    #("reads", None, None, "Reads to plot"),
    ("max_reads", 10, int, ""),
    ("show_events", False, bool, "Display event means plotted over signal instead of model current (currently events are computed from scratch)"),
    ("no_model", False, bool, "Will not plot the expected reference signal if True"),
    ("multi_background", False, bool, "Will plot multiple stacked background colors for multiple tracks if True"),
    ("yaxis_fixed", False, bool, ""),
    #("track_colors", ["purple", "darkgreen", "royalblue", "crimson"], list, ""),
    #("base_colors", ["#80ff80", "#6b93ff", "#ffbd00", "#ff8080"], list, "Colors for each base (A,C,G,T/U)"), 
    ("fill_layer", "base", str, ""),
    ("fill_track", 0, None, "")
)

class RefplotParams(config.ParamGroup):
    _name = "refplot"
RefplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("layer", "current", str, ""),
    ("plot", None, str, "Type of plot to display (hist, violin, box)"),
    ("kmer_coord", None, int, "Center position of k-mer to display"),
)

class TrackplotParams(config.ParamGroup):
    _name = "trackplot"
TrackplotParams._def_params(
    ("tracks", None, None, "DTW aligment tracks"),
    ("panels", None, None, "List of tuples specifying which panels to display. First element of each tuple specifies plot type (mat, box, line, scatter), second element specifies which layer(s) to display."),
    #("track_colors", ["#AA0DFE", "#1CA71C", "#4676FF", "red"], list, ""),
    ("select_ref", None, str, "Reference Selection"),
    ("select_read", None, str, "Read Selection"),
    ("hover_read", False, bool, "If True will display read_id in mat hover"),
    ("show_legend", True, bool, "If True will display legend"),
    ("share_reads", False, bool, "If True will only display reads shared by all alignment tracks with shared y-axis"),
    ("width", None, int, "Figure width"),
    ("panel_heights", None, None, "Relative height of each panel"),
    ("min_height", 200, int, "Minimum figure height"),
    ("outfile", None, str, "Output file"),
)

class ReadstatsParams(config.ParamGroup):
    _name = "readstats"
ReadstatsParams._def_params(
    ("stats", ["model_diff"], None, "Which statistics to compute and output"),
    ("pca_layer", "current", str, "Which statistics to use for PCA"),
    ("pca_components", 2, int, "Number of principle components to output for the \"pca\" command."),
    ("summary_stats", ["mean"], None, "Summary statistics to compute for \"model_diff\" command."),
)

