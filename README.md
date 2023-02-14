# Uncalled4

A **U**tility for **N**anopore **C**urrent **Al**ignment to **L**arge **E**xpanses of **D**NA

![UNCALLED logo](logo.png "UNCALLED logo")

A toolkit for nanopore signal alignment, analysis, and visualization

Features an alignment algorithm guided by Guppy metadata, methods for comparing Tombo and Nanopolish alignments, epigenetic modification detection statistics, and interactive alignment visualizations

Uncalled4 is under active development, and things like command line arguments and file formats may change in future versions. It has primarily been developed on r9.4 RNA data, and certain functionalities may not work for DNA yet.

For [real-time targeted sequencing](https://www.nature.com/articles/s41587-020-0731), see the [main branch](https://github.com/skovaka/UNCALLED)

## Table of Contents

- [Installation](#installation)
- [`index`: Reference Indexing](#index)
- [DTW Alignment and Storage](#dtw-alignment-and-storage)
  - [`dtw`: Perform DTW alignment guided by basecalled alignments](#dtw)
  - [`train`: Train new k-mer pore models](#train)
  - [`convert`: Import DTW alignments produced by Nanopolish or Tombo](#convert)
  - [`db`: Edit, merge, and ls alignment databases](#db)
- [DTW Visualization](#dtw-visualization)
  - [`dotplot`: Plot signal-to-reference alignment dotplots](#dotplot)
  - [`trackplot`: Plot alignment tracks and per-reference statistics](#trackplot)
  - [`browser`: Interactive signal alignment genome browser](#browser)
- [DTW Analysis](#dtw-analysis)
  - [`refstats`: Calculate per-reference-coordinate statistics](#refstats)
  - [`layerstats`: Compute, compare, and query alignment layers](#layerstats)
- [Alignment Layers](#alignment-layers)
- [Future Work](#future-work)
- [Release Notes](#release-notes)

## Installation

```
> pip install git+https://github.com/skovaka/UNCALLED.git@uncalled4
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

Build an index from a FASTA reference. Must be performed before any other analysis.

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

`--no-bwt` should be used if you don't plan to use this index for `uncalled realtime`, `map`, or `sim`

Note that UNCALLED uses the [BWA](https://github.com/lh3/bwa) FM Index to encode the reference, and this command will use a previously built BWA index if all the required files exist with the specified prefix. Otherwise, a new BWA index will be automatically built.

## DTW Alignment and Storage

Commands for generating, reading, and converting nanopore signal alignments are described in the following sections. Several input output formats are supported, specified by the flags `--<fmt>-in/--<fmt>-out`. Only one input and output can be specified per command, and not all inputs are supported for each command. Output formats are summarized below:

- `--sql-in/--sql-out`: Sqlite3 database format. This is the most versitile format, but also the slowest. It is currently the only format supported by all DTW-related subcommands, since it supports random access and storage of many datatypes. Multiple tracks can be stored in one database, and track names can be specified in the format `<file.db>:<track_name1>[,<track_name2>...]`
- `--bam-in/--bam-out`: Experimental storage of alignments in BAM tags. This is distinct from nanopolish's "--sam" option in that it stores all information in tags, preserving basecalled alignment information, and includes per-reference current levels so FAST5 files don't need to be accessed for most operations. Sorted and indexed BAM files can be loaded by `browser`, and we plan to eventually support BAM input for all commands.
- `--eventalign-in/--eventalign-out`: Nanopolish "eventalign" format. Can be used to import nanopolish alignments, or to use Uncalled4/Tombo alignments with tools which accept Nanopolish alignments. For `--eventalign-in`, read IDs (`--print-read-names`) and sample coordaintes (`signal-index`) must be included. Use the `--eventalign-flags` option to include these and other optional fields with `--eventalign-out`
- `--tsv-out`: General tab-seperated-values output. Can specify which layers to output with the "--tsv-cols" option (see [Alignment Layers](#alignment-layers) for options)
- `--tombo-in`: Currently only availible with the `convert` command to import Tombo alignments. This is the least tested input format, and currently only works for r9.4 RNA.

Note that only SQL and BAM output store configuation information within the file. For other formats, parameters like reference index, fast5 paths/index, and RNA status will need to specified when reading the file.

Below are summaries for each subcommand. See `uncalled <subcommand> -h` for more detailed documentation

### `dtw`

Perform DTW alignment guided by basecalled alignments

Currently the fast5 files must contain basecalling information output by Guppy via the `--fast5_out` option, or basecaller moves should be in the BAM input file with the Guppy `--moves_out` option.


```
usage: uncalled dtw [-h] --bam-in [BAM_IN] [-r] [-l READ_FILTER] [-x FAST5_INDEX] [-n MAX_READS]
                    [--sql-out SQL_OUT | --tsv-out [TSV_OUT] | --bam-out [BAM_OUT] |
                    --eventalign-out [EVENTALIGN_OUT]] [--tsv-cols TSV_COLS] [--tsv-na [TSV_NA]]
                    [--eventalign-flags EVENTALIGN_FLAGS] [--mask-skips [MASK_SKIPS]]  
                    [--mask-indels MASK_INDELS] [-f] [-a] [--bc-cmp] [-p PORE_MODEL]             
                    [--save-bands] [--full-overlap] [--rna] [-R REF_BOUNDS] [-i ITERATIONS]
                    [-c {abs_diff,z_score,norm_pdf}] [--skip-cost SKIP_COST]                     
                    [--stay-cost STAY_COST] [--move-cost MOVE_COST] [-b BAND_WIDTH]        
                    [-s BAND_SHIFT] [-N {ref_mom,model_mom}] [--bc-group BC_GROUP]               
                    index_prefix fast5_files [fast5_files ...]                  
```

### `train`

Iteratively train a new k-mer pore model

Accepts most of the same paramters as `uncalled dtw`, in addition to number of iterations. First iteration must use some starting pore model, while subsequent iterations use  the pore model from the previous iteration.

```
usage: uncalled train [-h] [-i ITERATIONS] [--kmer-samples KMER_SAMPLES] [--buffer-size BUFFER_SIZE]          
                      [-d MAX_BCALN_DIST] [--use-median] [--out-dir OUT_DIR] [-a] [--skip-dtw] [-p PROCESSES] 
                      [--bam-chunksize BAM_CHUNKSIZE] [--guppy-in GUPPY_IN] --bam-in [BAM_IN]                 
                      [--out-name OUT_NAME] [-r] [-l READ_FILTER] [-x FAST5_INDEX] [-n MAX_READS]             
                      [--del-max DEL_MAX] [--mask-skips [MASK_SKIPS]] [--mask-indels MASK_INDELS] [-f]        
                      [-m PORE_MODEL] [--kmer-shift KMER_SHIFT] [--save-bands] [--full-overlap] [--rna]       
                      [-R REF_BOUNDS] [-c {abs_diff,z_score,norm_pdf}] [--skip-cost SKIP_COST]                
                      [--stay-cost STAY_COST] [--move-cost MOVE_COST] [-b BAND_WIDTH] [-s BAND_SHIFT]         
                      [-N {ref_mom,model_mom}] [--norm-median] [--norm-seg] [--bc-group BC_GROUP] [-C CONFIG] 
                      index_prefix fast5_files [fast5_files ...]                                              
```

### convert

Convert between signal alignment file formats

```
usage: uncalled convert [-h] [--sql-in SQL_IN | --bam-in BAM_IN | --eventalign-in [EVENTALIGN_IN]
                        | --tombo-in TOMBO_IN]                                                   
                        [--sql-out SQL_OUT | --eventalign-out [EVENTALIGN_OUT] | --tsv-out
                        [TSV_OUT]] [--tsv-cols TSV_COLS] [--eventalign-flags EVENTALIGN_FLAGS]
                        [--mask-skips [MASK_SKIPS]] [--fast5s FAST5S [FAST5S ...]]          
                        [-l READ_FILTER] [-x FAST5_INDEX] [-r] [--rna] [-R REF_BOUNDS] [-f] [-a]
                        index_prefix       
```
                                                                                              

### `db`

Edit, merge, and ls SQL alignment databases. Only applicable to files generated via `--sql-out`.

```
uncalled db [-h]

subcommand options:
ls       List all tracks in a database
delete   Delete a track from a database
merge    Merge databases into a single file
edit     Rename, change fast5 paths, or set description
```

## DTW Visualization

All visualizations are generated using Plotly.

### `dotplot`

Plot signal-to-reference alignment dotplots

If any alignments were generated using `uncalled dtw`, will also display a dotplot of the basecalled alignment projected into signal space in orange.

```
usage: uncalled dotplot [-h] [--sql-in SQL_IN | --bam-in [BAM_IN] | --eventalign-in
                        [EVENTALIGN_IN]] [-o OUT_PREFIX] [--ref REF]                   
                        [--fast5s FAST5S [FAST5S ...]] [-x FAST5_INDEX] [-r] [--rna]
                        [-f {png,pdf,svg}] [-R REF_BOUNDS] [-l READ_FILTER] [-L LAYERS]
                        [-b BCALN_TRACK] [--multi-background] [--show-events] [--show-bands]
                        [--no-model] [--bcaln-error]                                         
```

### `trackplot`

Plot alignment tracks and per-reference statistics

Trackplots are defined by a series of panels displaying different layers. A `mat` panel display a heatmap of layer values for each ref/read coordinate on each track. A `box` panel displays boxplots of layer summary statistics for each track. `line` and `scatter` panels display [`refstats`](#refstats) summary statistics, specified by `<layer>.<statistic>` (e.g. `current.mean`, `model_diff.median`).

```
usage: uncalled trackplot [-h] [--sql-in SQL_IN | --bam-in [BAM_IN]] [--ref REF]         
                          [--fast5s FAST5S [FAST5S ...]] [-x FAST5_INDEX] [-r] [--rna]
                          [--pore-model PORE_MODEL] [-f] [-l READ_FILTER]                  
                          [-H PANEL_HEIGHTS [PANEL_HEIGHTS ...]] [--shared-refs-only]
                          [--shared-reads-only] [--share-reads] [--hover-read] [-o OUTFILE]     
                          [--mat LAYER] [--box LAYER] [--line LAYER.STAT] [--scatter LAYER.STAT]
                          ref_bounds           
```

### `browser`

Interactive signal alignment genome browser

This feature is in very early stages. Currently it features trackplot visualization where you can click on different locations to display read/reference information.

```
usage: uncalled browser [-h] [--sql-in SQL_IN | --bam-in [BAM_IN]] [--ref REF]
                        [--fast5s FAST5S [FAST5S ...]] [-x FAST5_INDEX] [-r] [--rna]
                        [-l READ_FILTER] [-f] [--pore-model PORE_MODEL] [-p BROWSER_PORT]
                        [-o OUTFILE]
                        ref_bounds
```

## DTW Analysis

These functions compute statistics over reference and read coordinates. `refstats` computes summary and comparison statistics (e.g. Kolmogorov-Smirnov test) over reference coordinates, while `layerstats` maintains the read-by-reference dimensions of the DTW alignments.

(Coming soon: `readstats` to compute read-level statistics)

### `refstats`

Calculate per-reference-coordinate statistics

```
uncalled refstats [-R REF_BOUNDS] [-C REF_CHUNKSIZE] [-c] [-v]
                  input [input ...] refstats_layers refstats

positional arguments:
  input                 Input tracks specifier. Should be in the format
                        <file.db>[:<track1>[,<track2>...]]. If no track
                        names are specified, all tracks will be loaded   from
                        the database.

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
```

### `layerstats`

This subcommand features `compare`, which measures distance between alignments of the same set of reads (e.i. different alignment methods). It also features `dump`, which outputs DTW layers from a database.

#### `compare`

Compute distance between alignments of the same reads

Note: currently the `--save` option will only work if no other comparisons have been saved to the primary track (the first specified track) 

```
uncalled layerstats compare [-b] [-s] input [input ...]

positional arguments:
input               Input tracks specifier. Should be in the format
                    <file.db>[:<track1>[,<track2>...]]. If no track
                    names are specified, all tracks will be loaded from
                    the database.


optional arguments:
  -b, --bcaln  Compare against basecalled alignment. If two tracks input will
               look for "bcaln" group in second track, otherwise will look in
               the first track. (default: False)
  -s, --save   Will save in database if included, otherwise will output
               to TSV. Will be associated with the first track listed.
               (default: False)
```

## Alignment Layers

Uncalled4 stores signal alignments as a set of **layers** associated with read and reference coordinates. Each layer is derived from the read signal (e.g. `current`), the reference sequence (e.g. `kmer`), or other layers (e.g. `model_diff`). Layers are organized into **layer groups**: `dtw` for signal alignments, `bcaln` for projected basecalled alignments, and `cmp` for alignment comparisons. Several subcommands detailed above take layer names as input, which should generaly be in the form `<group>.<layer>`. For brevity, group can be excluded for `dtw` layers (e.g. you can simply input `current`, not `dtw.currnt`). Below is a table of currently available layers:

| Group | Layer   | Description |
|-------|---------|-------------|
| dtw   | current | Mean read signal current (pA) |
| dtw   | start | Read signal sample start index |
| dtw   | length | Read signal sample length |
| dtw   | dwell   | Signal dwell time (ms/nt, proportional to length) |
| dtw   | model_diff | Difference between predicted (via a pore model) and observed current|
| dtw   | abs_diff | Absolute value of dtw.model_diff |
| dtw   | kmer | Binarized reference k-mer |
| dtw   | events | Number of raw signal events aligned ("stays" > 1, "skips" < 1) |
| dtw   | base | Binarized reference base |
| bcaln   | start | Estimated read signal sample start index |
| bcaln   | length | Estimated read signal sample length |
| bcaln   | indel | Number of inserted (>0) or deleted (<0) nucleotides |
| cmp   | mean_ref_dist | Mean reference distance between two alignments (must first run [`layerstats compare`](#compare)) |
| cmp   | jaccard | Raw sample jaccard distances between two alignments (must first run [`layerstats compare`](#compare)) |

Read alignments are grouped into **tracks**, which can be stored in an sqlite3 database, or individual BAM, eventalign, or TSV files. When running [`uncalled dtw`](#dtw) or [`uncalled convert`](#convert) with the `--sql-out` option, you can specify an output database file along with an optional track name in the format `<filename.db>:<track_name>`. If only `<filename.db>` is specified, the track name will be the filename without the extension. Multiple tracks can be input as comma-separated track names: `<filename.db>:<track1>,<track2>,...`. You can change the name of a track with [`uncalled db edit`](#db), and also set a human readable description to appear in visualizations.

All tracks must be written to the same database for multi-track visualization and analysis (e.g. comparing alignments, calculating KS statistics). You can merge multiple databases into a single file using [`uncalled db merge`](#db)

## Future Work

Uncalled4 (v3.3.0) is a work in progress. Many additional features and optimizations are planned, which may require changes to the command line interface or database format. We also plan to provide a Python API in the future.

Real-time targeted sequencing (`uncalled realtime`) is currently unavailable in Uncalled4. Related subcommands (`map`, `sim`, etc.) are available, but aren't documented here for brevity. We plan to integrate all functionalities eventually.

## Release Notes
- v3.3: added BAM input/output, generalized TSV output, and started accepting non-SQL input (BAM, eventalign) for some commands. Changed all alignment input/outputs to optional arguments (not positional), and refactored `convert` command structure. DTW now must be guided by BAM files **no longer supporting PAF files for `dtw`**
- v3.2: expanded k-mer model to begin supporting R10 alignment
- v3.1: introduced Plotly visualizations and sqlite3 database
- v3.0: added DTW alignment, analysis, and visualization commands
- v2.2: added event profiler which masks out pore stalls, and added compile-time debug options
- v2.1: updated ReadUntil client for latest MinKNOW version, made `uncalled index` automatically build the BWA index, added hdf5 submodule, further automated installation by auto-building hdf5, switched to using setuptools, moved submodules to submods/
- v2.0: released the ReadUntil simulator `uncalled sim`, which can predict how much enrichment UNCALLED could provide on a given reference, using a control and UNCALLED run as a template. _Also changed the format of certain arguments_: index prefix and fast5 list are now positional, and some flags have changed names.
- v1.2: fixed indexing for particularly large or small reference
- v1.1: added support for altering chunk size
- v1.0: pre-print release
