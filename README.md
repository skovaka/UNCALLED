# Uncalled4

A **U**tility for **N**anopore **C**urrent **Al**ignment to **L**arge **E**xpanses of **D**NA

![UNCALLED logo](uncalled_logo_small.png "UNCALLED logo")

A toolkit for nanopore signal alignment, analysis, and visualization

Features an alignment algorithm guided by Guppy metadata, methods for comparing Tombo and Nanopolish alignments,  epigenetic modification detection statistics, and interactive alignment visualizations

For [real-time targeted sequencing](https://www.nature.com/articles/s41587-020-0731), see the [main branch](https://github.com/skovaka/UNCALLED)

## Table of Contents

- [Installation](#installation)
- [`index`: Reference Indexing](#index)
- [DTW Alignment and Storage](#dtw-alignment-and-storage)
 - [`dtw`: Perform DTW alignment guided by basecalled alignments](#dtw)
 - [`convert`: Import DTW alignments produced by Nanopolish or Tombo](#convert)
 - [`db`: Edit and query alignment database metadata](#db)
- [DTW Visualization](#dtw-visualization)
 - [`dotplot`: Plot signal-to-reference alignment dotplots](#dotplot)
 - [`trackplot`: Plot alignment tracks and per-reference statistics](#trackplot)
 - [`browser`: Interactive signal alignment genome browser](#browser)
- [DTW Analysis](#dtw-analysis)
 - [`refstats`: Calculate per-reference-coordinate statistics](#refstats)
 - [`dtwstats`: Compute, compare, and query alignment layers](#dtwstats)
- [Release Notes](#release-notes)

## Installation

```
> pip install git+https://github.com/skovaka/UNCALLED.git@uncalled4 --user
```

OR

```
> git clone --recursive -b uncalled4 https://github.com/skovaka/UNCALLED.git
> cd UNCALLED
> pip install .
```

Requires python >= 3.8

Other dependencies are included via submodules, so be sure to clone with `git --recursive`

Uncalled4 is compatible with Linux and Mac OSX. It has been primarily developed and tested on Ubuntu Linux.

## `index`

Build an index from a FASTA reference

```
uncalled index [-o INDEX_PREFIX] [--no-bwt] fasta_filename

positional arguments:
  fasta_filename        FASTA file to index

optional arguments:
  -o INDEX_PREFIX, --index-prefix INDEX_PREFIX
                        Index output prefix. Will use input fasta filename by
                        default (default: None)
  --no-bwt              Will only generate the pacseq if specified, which is
                        much faster to build. Can only be used with DTW
                        subcommands (NOT map, sim, or realtime) (default:
                        False)
```

Note that UNCALLED uses the [BWA](https://github.com/lh3/bwa) FM Index to encode the reference, and this command will use a previously built BWA index if all the required files exist with the specified prefix. Otherwise, a new BWA index will be automatically built.

## DTW Alignment and Storage

The following subcommands generate and update dynamic time warping (DTW) alignment tracks. Alignment tracks are currently stored in a sqlite3 database. Multiple tracks can and should be stored in the same database in order to be analyzed together. Currently there is no way to merge databases. Tracks can be specified in `<file.db>:<track_name>` format.

### `dtw`

Perform DTW alignment guided by basecalled alignments

Currently the **fast5 files must contain basecalling information** output by Guppy via the `--fast5_out` option.


```
uncalled dtw [-r] [-l READ_FILTER] [-n MAX_READS]
             -m MM2_PAF [-o OUTPUT][--full-overlap]
              [--rna] [-R REF_BOUNDS] [-a] [-f]
              index_prefix fast5_files [fast5_files ...]


positional arguments:
  index_prefix          BWA index prefix
  fast5_files           List of paths to any combination of: 1. fast5 files 2.
                        directories to search for fast5 files (optionally
                        recursive) 3. text file(s) listing one fast5 file or
                        directory to search per line

optional arguments:
  --rna                 Should be set for direct RNA data
  -m MM2_PAF, --mm2-paf MM2_PAF
                        Path to minimap2 alignments of basecalled reads in PAF
                        format. Used to determine where each should be
                        aligned. Should include cigar string. (default: None)
  -r, --recursive       Will search directories recursively for any file                       
                        ending in ".fast5" if true (default: False)
  -l READ_FILTER, --read-filter READ_FILTER
                        List of read IDs and/or text file(s) containing one
                        read ID per line (default: [])
  -n MAX_READS, --max-reads MAX_READS
                        Maximum number of reads to load. (default: 0)
  -o OUTPUT, --output OUTPUT
                        Output track specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is the output
                        sqlite database. If <track_name> is not specified,
                        will be same as filename (without extension) (default:
                        None)
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates
                        (default: None)
  --full-overlap        If true will only include reads which fully cover
                        reference bounds (default: False)
  -f, --overwrite       Overwrite existing tracks (default: False)
  -a, --append          Append reads to existing tracks (default: False)
```


The `uncalled dtw` command consists of several subcommands for alignment, analysis, and visualizations of raw nanopore signal aligned via dynamic time warping. `uncalled dtw align` and `uncalled dtw convert` produce raw signal **alignment tracks**, which is a set of reads (typically from one sample) aligned to the same reference. The rest of the subcommands detailed below produce visualizations, analysis, and comparisons of alignment tracks.

### convert

Import DTW alignments produced by Nanopolish and Tombo

#### `nanopolish`

Currently nanopolish eventalign **must** but run with the options `-n/--print-read-names`, `--signal-index`, and `--scale-events`

```
uncalled nanopolish [-r] [-l READ_FILTER] [-n MAX_READS] [--rna]
                           [-R REF_BOUNDS] [-f] [-a] -o OUTPUT -x FAST5_INDEX
                           index_prefix fast5_files [fast5_files ...]
                           eventalign_tsv

positional arguments:
  index_prefix          BWA index prefix
  fast5_files           List of paths to any combination of: 1. fast5 files 2.
                        directories to search for fast5 files (optionally
                        recursive) 3. text file(s) listing one fast5 file or
                        directory to search per line
  eventalign_tsv

optional arguments:

  -x FAST5_INDEX, --fast5-index FAST5_INDEX
                        Filename mapping between read IDs and fast5 files Can
                        be sequencing summary output by basecaller,
                        "filename_mapping.txt" from ont_fast5_api, or
                        nanopolish "*.index.readdb" file (default: )
  -r, --recursive       Will search directories recursively for any file
                        ending in ".fast5" if true (default: False)
  -l READ_FILTER, --read-filter READ_FILTER
                        List of read IDs and/or text file(s) containing one
                        read ID per line (default: [])
  -n MAX_READS, --max-reads MAX_READS
                        Maximum number of reads to load. (default: 0)
  --rna                 Should be set for direct RNA data
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates
                        (default: None)
  -f, --overwrite       Overwrite existing tracks (default: False)
  -a, --append          Append reads to existing tracks (default: False)
  -o OUTPUT, --output OUTPUT
                        Output track specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is the output
                        sqlite database. If <track_name> is not specified,
                        will be same as filename (without extension) (default:
                        None)
```

#### `tombo`

```
uncalled tombo [-h] [-r] [-l READ_FILTER] [-n MAX_READS] [--rna]
                      [-R REF_BOUNDS] [-f] [-a] -o OUTPUT
                      index_prefix fast5_files [fast5_files ...]

Convert from Tombo resquiggled fast5s to uncalled DTW track

positional arguments:
  index_prefix          BWA index prefix
  fast5_files           List of paths to any combination of: 1. fast5 files 2.
                        directories to search for fast5 files (optionally
                        recursive) 3. text file(s) listing one fast5 file or
                        directory to search per line

optional arguments:
  -h, --help            show this help message and exit
  -r, --recursive       Will search directories recursively for any file
                        ending in ".fast5" if true (default: False)
  -l READ_FILTER, --read-filter READ_FILTER
                        List of read IDs and/or text file(s) containing one
                        read ID per line (default: [])
  -n MAX_READS, --max-reads MAX_READS
                        Maximum number of reads to load. (default: 0)
  --rna                 Should be set for direct RNA data
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates
                        (default: None)
  -f, --overwrite       Overwrite existing tracks (default: False)
  -a, --append          Append reads to existing tracks (default: False)
  -o OUTPUT, --output OUTPUT
                        Output track specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is the output
                        sqlite database. If <track_name> is not specified,
                        will be same as filename (without extension) (default:
                        None)
```

### `db`

Edit and query alignment database metadata

```
usage: uncalled db [-h]

subcommand options:
ls       List all tracks in a database
delete   Delete a track from a database
edit     Rename, change fast5 paths, or set description
```

## DTW Visualization

All visualizations are generated using Plotly.

### `dotplot`

Plot signal-to-reference alignment dotplots

If any alignments were generated using `uncalled dtw`, will also display a dotplot of the basecalled alignment projected into signal space in orange.

```
usage: uncalled dotplot [-h] [-o OUT_PREFIX] [-f {png,svg,pdf}]
                        [-R REF_BOUNDS] [-l READ_FILTER] [-L LAYERS]
                        input [input ...]

positional arguments:
  input                 Input tracks specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is an
                        Uncalled4 aligment track database and <track_name>
                        optionally specifies which tracks to read (reads all
                        by default)

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PREFIX, --out-prefix OUT_PREFIX
                        If included will output images with specified prefix,
                        otherwise will display interactive plot. (default:
                        None)
  -f {png,svg,pdf}, --out-format {png,svg,pdf}
                        Image output format. Only has an effect with -o
                        option. (default: svg)
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates
                        (default: None)
  -l READ_FILTER, --read-filter READ_FILTER
                        Only load reads which overlap these coordinates
                        (default: None)
  -L LAYERS, --layers LAYERS
```

### `trackplot`

Plot alignment tracks and per-reference statistics

Trackplots are defined by a series of panels displaying different layers. A `mat` panel display a heatmap of layer values for each ref/read coordinate on each track. A `box` panel displays boxplots of layer summary statistics for each track. `line` and `scatter` panels display [`refstats`](#refstats) summary statistics, specified by `<layer>.<statistic>`.

```
usage: uncalled trackplot [-h] [-f] [-H PANEL_HEIGHTS [PANEL_HEIGHTS ...]]
                          [--mat LAYER] [--box LAYER] [--line LAYER.STAT]
                          [--scatter LAYER.STAT] [-o OUTFILE]
                          input [input ...] ref_bounds


positional arguments:
  input                 Input tracks specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is an
                        Uncalled4 aligment track database and <track_name>
                        optionally specifies which tracks to read (reads all
                        by default)
  ref_bounds            Only load reads which overlap these coordinates

optional arguments:
  -h, --help            show this help message and exit
  -f, --full-overlap    If true will only include reads which fully cover
                        reference bounds (default: False)
  -H PANEL_HEIGHTS [PANEL_HEIGHTS ...], --panel-heights PANEL_HEIGHTS [PANEL_HEIGHTS ...]
                        Relative height of each panel (default: None)
  --mat LAYER           Display a ref-by-read matrix of specified alignment
                        layer (default: None)
  --box LAYER           Display a boxplot of specified layer (default: None)
  --line LAYER.STAT     Display a line plot of specifed layer summary
                        statistic (default: None)
  --scatter LAYER.STAT  Display a line plot of specifed layer summary
                        statistic (default: None)
  -o OUTFILE, --outfile OUTFILE
                        Output file (default: None)
```

### `browser`

Interactive signal alignment genome browser

This feature is in very early stages. Currently it features trackplot visualization where you can click on different locations to display read/reference information.

```
uncalled browser [-r REFSTATS] [-f] [-o OUTFILE]
                 input [input ...] ref_bounds


positional arguments:
  input                 Input tracks specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is an
                        Uncalled4 aligment track database and <track_name>
                        optionally specifies which tracks to read (reads all
                        by default)
  ref_bounds            Only load reads which overlap these coordinates

optional arguments:
  -r REFSTATS, --refstats REFSTATS
                        Per-reference summary statistics to compute for each
                        layer (default: None)
  -f, --full-overlap    If true will only include reads which fully cover
                        reference bounds (default: False)
  -o OUTFILE, --outfile OUTFILE
                        Output file (default: None)
```

## DTW Analysis

These functions compute statistics over reference and read coordinates. `refstats` computes summary and comparison statistics (e.g. Kolmogorov-Smirnov test) over reference coordinates, while `dtwstats` maintains the read-by-reference dimensions of the DTW alignments.

(Coming soon: `readstats` to compute read-level statistics)

### `refstats`

Calculate per-reference-coordinate statistics

```
uncalled refstats [-R REF_BOUNDS] [-C REF_CHUNKSIZE] [-c] [-v]
                  input [input ...] refstats_layers refstats

positional arguments:
  input                 Input tracks specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is an
                        Uncalled4 aligment track database and <track_name>
                        optionally specifies which tracks to read (reads all
                        by default)
  refstats_layers       Comma-separated list of layers over which to compute
                        summary statistics
  refstats              Comma-separated list of summary statistics to compute.
                        Some statisitcs (ks) can only be used if exactly two
                        tracks are provided {ks,mean,q95,median,skew,max,var,m
                        in,q5,stdv,kurt,q25,q75}
optional arguments:
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates
                        (default: None)
  -C REF_CHUNKSIZE, --ref-chunksize REF_CHUNKSIZE
                        Number of reference coordinates to query for iteration
                        (default: 10000)
  -c, --cov             Output track coverage for each reference position (default: False)
  -v, --verbose-refs    Output reference name and strand (default: False)
```

### `dtwstats`

This subcommand features `compare`, which measures distance between alignments of the same set of reads (e.i. different alignment methods). It also features `dump`, which outputs DTW layers from a database.

#### `compare`

Compute distance between alignments of the same reads

```
uncalled dtwstats compare [-b] [-s] input [input ...]

positional arguments:
  input        Input tracks specifier. Should be in format
               <file.db>[:<track_name>], where file.db is an Uncalled4
               aligment track database and <track_name> optionally specifies
               which tracks to read (reads all by default)

optional arguments:
  -b, --bcaln  Compare against basecalled alignment. If two tracks input will
               look for "bcaln" group in second track, otherwise will look in
               the first track. (default: False)
  -s, --save   Will save in database if included, otherwise will output to TSV
               (default: False)
```

#### `dump`

Output DTW alignment paths and statistics

```
uncalled dump [-R REF_BOUNDS] [-l READ_FILTER] [-o {db,tsv}]
              input [input ...] layers [layers ...]

positional arguments:
  input                 Input tracks specifier. Should be in format
                        <file.db>[:<track_name>], where file.db is an
                        Uncalled4 aligment track database and <track_name>
                        optionally specifies which tracks to read (reads all
                        by default)
  layers                Layers to retrieve or compute

optional arguments:
  -h, --help            show this help message and exit
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates
                        (default: None)
  -l READ_FILTER, --read-filter READ_FILTER
                        Only load reads which overlap these coordinates
                        (default: None)
```


## Release notes
- v3.1: introduced Plotly visualizations and sqlite3 database
- v3.0: added DTW alignment, analysis, and visualization commands
- v2.2: added event profiler which masks out pore stalls, and added compile-time debug options
- v2.1: updated ReadUntil client for latest MinKNOW version, made `uncalled index` automatically build the BWA index, added hdf5 submodule, further automated installation by auto-building hdf5, switched to using setuptools, moved submodules to submods/
- v2.0: released the ReadUntil simulator `uncalled sim`, which can predict how much enrichment UNCALLED could provide on a given reference, using a control and UNCALLED run as a template. _Also changed the format of certain arguments_: index prefix and fast5 list are now positional, and some flags have changed names.
- v1.2: fixed indexing for particularly large or small reference
- v1.1: added support for altering chunk size
- v1.0: pre-print release
