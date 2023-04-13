# Uncalled4

A **U**tility for **N**anopore **C**urrent **Al**ignment to **L**arge **E**xpanses of **D**NA

![UNCALLED logo](logo.png "UNCALLED logo")

A toolkit for nanopore signal alignment, analysis, and visualization

Features an alignment algorithm guided by Guppy metadata, methods for comparing Tombo and Nanopolish alignments, epigenetic modification detection statistics, and interactive alignment visualizations

Uncalled4 is under active development, and things like command line arguments and file formats may change in future versions.

For [real-time targeted sequencing](https://www.nature.com/articles/s41587-020-0731), see the [main branch](https://github.com/skovaka/UNCALLED)

## Table of Contents

- [Installation](#installation)
- [DTW Alignment and Storage](#dtw-alignment-and-storage)
  - [`dtw`: Perform DTW alignment guided by basecalled alignments](#dtw)
  - [`train`: Train new k-mer pore models](#train)
  - [`convert`: Import DTW alignments produced by Nanopolish or Tombo](#convert)
- [DTW Visualization](#dtw-visualization)
  - [`trackplot`: Plot alignment tracks and per-reference statistics](#trackplot)
  - [`browser`: Interactive signal alignment genome browser](#browser)
- [DTW Analysis](#dtw-analysis)
  - [`refstats`: Calculate per-reference-coordinate statistics](#refstats)
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

## DTW Alignment and Storage

Commands for generating, reading, and converting nanopore signal alignments are described in the following sections. Several input output formats are supported, specified by the flags `--<fmt>-in/--<fmt>-out`. Only one input and output can be specified per command, and not all inputs are supported for each command. Output formats are summarized below:

- `--bam-in/--bam-out`: Uncalled4 primarily stores signal alignments in BAM alignment tags, alongside the standard basecalled read alignments. BAM files must be sorted and indexed prior to viewing within the browser. They can be converted to eventalign or TSV format via the convert command.
- `--eventalign-in/--eventalign-out`: Nanopolish "eventalign" format. Can be used to import nanopolish alignments, or to use Uncalled4/Tombo alignments with tools which accept Nanopolish alignments. For `--eventalign-in`, read IDs (`--print-read-names`) and sample coordaintes (`signal-index`) must be included. Use the `--eventalign-flags` option to include these and other optional fields with `--eventalign-out`
- `--tsv-out`: General tab-seperated-values output. Can specify which layers to output with the "--tsv-cols" option (see [Alignment Layers](#alignment-layers) for options)
- `--tombo-in`: Currently only availible with the `convert` command to import Tombo alignments. This is the least tested input format, and currently only works for r9.4 RNA.

Uncalled4 supports FAST5, SLOW5/BLOW5, and POD5 read signal formats. All `read_paths` will be searched for files ending with `.fast5`, `.slow5`, `.blow5`, or `.pod5`, optionally recursively. Uncalled4 will also attempt to automatically detect the appropriate pore model based on read file metadata, however this feature is experimental.

Below are summaries for each subcommand. See `uncalled <subcommand> -h` for more detailed documentation

### `dtw`

Perform DTW alignment guided by basecalled alignments

`ref_index` must be a FASTA file indexed via `samtools faidx`. `read_files` must contain at least one FAST5, SLOW5, BLOW5, or POD5 file, optionally recursively with the `--recursive` option.

Currently the fast5 files must contain basecalling information output by Guppy via the `--fast5_out` option, or basecaller moves should be in the BAM input file with the Guppy `--moves_out` option.


```
usage: uncalled dtw [-h] --bam-in [BAM_IN] [-r] [-l READ_FILTER] [-x READ_INDEX] [-n MAX_READS]
                    [--tsv-out [TSV_OUT] | --bam-out [BAM_OUT] |
                    --eventalign-out [EVENTALIGN_OUT]] [--tsv-cols TSV_COLS] [--tsv-na [TSV_NA]]
                    [--eventalign-flags EVENTALIGN_FLAGS] [--mask-skips [MASK_SKIPS]]  
                    [--mask-indels MASK_INDELS] [-f] [-a] [--bc-cmp] [-p PORE_MODEL]             
                    [--save-bands] [--full-overlap] [--rna] [-R REF_BOUNDS] [-i ITERATIONS]
                    [-c {abs_diff,z_score,norm_pdf}] [--skip-cost SKIP_COST]                     
                    [--stay-cost STAY_COST] [--move-cost MOVE_COST] [-b BAND_WIDTH]        
                    [-s BAND_SHIFT] [-N {ref_mom,model_mom}] [--bc-group BC_GROUP]               
                    ref_index read_files [read_files ...]                  
```

### `train`

Iteratively train a new k-mer pore model

Accepts most of the same paramters as `uncalled dtw`, in addition to number of iterations. First iteration must use some starting pore model, while subsequent iterations use  the pore model from the previous iteration.

```
usage: uncalled train [-h] [-i ITERATIONS] [--kmer-samples KMER_SAMPLES] [--buffer-size BUFFER_SIZE]          
                      [-d MAX_BCALN_DIST] [--use-median] [--out-dir OUT_DIR] [-a] [--skip-dtw] [-p PROCESSES] 
                      [--bam-chunksize BAM_CHUNKSIZE] [--guppy-in GUPPY_IN] --bam-in [BAM_IN]                 
                      [--out-name OUT_NAME] [-r] [-l READ_FILTER] [-x READ_INDEX] [-n MAX_READS]             
                      [--del-max DEL_MAX] [--mask-skips [MASK_SKIPS]] [--mask-indels MASK_INDELS] [-f]        
                      [-m PORE_MODEL] [--kmer-shift KMER_SHIFT] [--save-bands] [--full-overlap] [--rna]       
                      [-R REF_BOUNDS] [-c {abs_diff,z_score,norm_pdf}] [--skip-cost SKIP_COST]                
                      [--stay-cost STAY_COST] [--move-cost MOVE_COST] [-b BAND_WIDTH] [-s BAND_SHIFT]         
                      [-N {ref_mom,model_mom}] [--norm-median] [--norm-seg] [--bc-group BC_GROUP] [-C CONFIG] 
                      ref_index read_files [read_files ...]                                              
```

### convert

Convert between signal alignment file formats

```
usage: uncalled convert [-h] [--bam-in BAM_IN | --eventalign-in [EVENTALIGN_IN]
                        | --tombo-in TOMBO_IN]                                                   
                        [--eventalign-out [EVENTALIGN_OUT] | --tsv-out
                        [TSV_OUT]] [--tsv-cols TSV_COLS] [--eventalign-flags EVENTALIGN_FLAGS]
                        [--mask-skips [MASK_SKIPS]] [--reads READS [READS ...]]          
                        [-l READ_FILTER] [-x READ_INDEX] [-r] [--rna] [-R REF_BOUNDS] [-f] [-a]
                        ref_index       
```
                                                                                              
## DTW Visualization

All visualizations are generated using Plotly.

### `trackplot`

Plot alignment tracks and per-reference statistics

Trackplots are defined by a series of panels displaying different layers. A `mat` panel display a heatmap of layer values for each ref/read coordinate on each track. A `box` panel displays boxplots of layer summary statistics for each track. `line` and `scatter` panels display [`refstats`](#refstats) summary statistics, specified by `<layer>.<statistic>` (e.g. `current.mean`, `model_diff.median`).

```
usage: uncalled trackplot [-h] [--bam-in [BAM_IN]] [--ref REF]         
                          [--reads READS [READS ...]] [-x READ_INDEX] [-r] [--rna]
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
usage: uncalled browser [-h] --bam-in [BAM_IN] [--ref REF]
                        [--reads READS [READS ...]] [-x READ_INDEX] [-r] [--rna]
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

All tracks must be written to the same database for multi-track visualization and analysis (e.g. comparing alignments, calculating KS statistics). You can merge multiple databases into a single file using [`uncalled db merge`](#db)

## Future Work

Uncalled4 (v3.4.0) is a work in progress. Many additional features and optimizations are planned, which may require changes to the command line interface or database format. We also plan to provide a Python API in the future.

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
