# UNCALLED

A **U**tility for **N**anopore **C**urrent **Al**ignment to **L**arge **E**xpanses of **D**NA

![UNCALLED logo](uncalled_logo_small.png "UNCALLED logo")

A streaming algorithm for mapping raw nanopore signal to DNA references

Enables real-time enrichment or depletion on Oxford Nanopore Technologies (ONT) MinION runs via ReadUntil

Also supports standalone signal mapping of fast5 reads

## Installation

```
> git clone --recursive https://github.com/skovaka/UNCALLED.git
> cd UNCALLED
> python setup.py install
```

Most dependecies included via submodules, so be sure to clone with `git --recursive`

UNCALLED must be installed into a python environment. To install without root privileges use the `--user` or `--prefix=<local-dir>` flag when installing, or use a tool such as [virtualenv](https://virtualenv.pypa.io) or [anaconda](https://www.anaconda.com).

Requires python-dev and GCC >= 4.8.1

[HDF5](https://www.hdfgroup.org/downloads/hdf5/) must be installed. Libraries and headers should be in system paths (ie `$LD_LIBRARY_PATH` and `$CPATH` respectively), or specified by running `python setup.py build_ext --library-dirs <hdf5-location>/lib --include-dirs <hdf5-location>/include` prior to installation.

The use of ReadUntil currently requires an ONT developer licence. Contact [ONT support](https://nanoporetech.com/contact) for more information. UNCALLED must then be installed into the MinKNOW environment alongside the ReadUntil API, using the same command specified in the ReadUntil API private GitHub repository.

We recommend running on a Linux machine. UNCALLED has been successfully installed and run on Mac computers, but real-time ReadUntil has not been tested on a Mac. Installing UNCALLED has not been attempted on Windows.

## Indexing

**Example:**

```
> bwa index -p E.coli E.coli.fasta
> uncalled index -i E.coli.fasta -x E.coli
```

UNCALLED requires a [BWA](https://github.com/lh3/bwa) index. You can use a previously built BWA index, or build a new one with the BWA instance provided in the `bwa/` submodule.

Before aligning, certain reference-specific parameters must be computed using `uncalled index`. The `<fasta-reference>` should be the same FASTA file which was used to build the BWA index. This will create an additional file in the same directory as the BWA index named `<bwa-prefix>.uncl`.

## Fast5 Mapping

**Example:**

```
> uncalled map -t 16 -x E.coli -i fast5_list.txt > uncalled_out.paf
Loading fast5s
Mapping

> head -n 4 uncalled_out.paf
b84a48f0-9e86-47ef-9d20-38a0bded478e 3735 77 328 + Escherichia_coli_chromosome 4765434 2024611 2024838 66 228 255  ch:i:427 st:i:50085  mt:f:53.662560
77fe7f8c-32d6-4789-9d62-41ff482cf890 5500 94 130 + Escherichia_coli_chromosome 4765434 2333754 2333792 38 39  255  ch:i:131 st:i:238518 mt:f:19.497091
eee4b762-25dd-4d4a-8a59-be47065029be 2905     *      *      *      *      *      *      *      *      *       255  ch:i:44  st:i:302369 mt:f:542.985229
e175c87b-a426-4a3f-8dc1-8e7ab5fdd30d 8052 84 154 + Escherichia_coli_chromosome 4765434 1064550 1064614 41 65  255  ch:i:182 st:i:452368 mt:f:38.611683
```

Arguments:

- `-x/--bwa-prefix` the prefix of the index to align to. Should be a BWA index that `uncalled index` was run on
- `-i/--fast5-files`  a text file containing the path to one fast5 file per line
- `-t/--threads` number of threads to use for mapping (default: 1)
- `-n/--read-count` maximum number of reads to map
- `-f/--filter` text file containing subset of read IDs (one per line) to map from the fast5 files (will map all by default)
- `-e/--max-events-proc` number of events to attempt mapping before giving up on a read (default 30,000). Note that there are approximately two events per nucleotide on average.


See [example/](example/) for a simple read and reference example.

## Real-Time ReadUntil

**Example:**

```
> uncalled list-ports
MN02686 (2019-11-18 12:30:56): 8000

> /opt/ONT/MinKNOW/ont-python/bin/uncalled realtime --port 8000 -t 16 -x E.coli --enrich -c 3 --post-script basecall.sh > uncalled_out.paf 
Starting client
Starting mappers
Mapping
Running post-script: 'basecall.sh'

> head -n 4 uncalled_out.paf
81ba344d-60df-4688-b37f-9064e76a3eb8 1352 *     *     *     *      *      *      *      *      *   255 ch:i:68  st:i:29101 mt:f:375.93841 wt:f:1440.934 mx:f:0.152565
404113c1-6ace-4690-885c-9c4a47da6476 450  *     *     *     *      *      *      *      *      *   255 ch:i:106 st:i:29268 mt:f:63.272270 wt:f:1591.070 en:f:0.010086
d9acafe3-23dd-4a0f-83db-efe299ee59a4 1355 *     *     *     *      *      *      *      *      *   255 ch:i:118 st:i:29378 mt:f:239.50201 wt:f:1403.641 ej:f:0.120715
8a6ec472-a289-4c50-9a75-589d7c21ef99 451  98 369 + Escherichia_coli 4765434 3421845 3422097 56 253 255 ch:i:490 st:i:29456 mt:f:79.419411 wt:f:8.551202 kp:f:0.097424
```

As mentioned above, you must sign an ONT developer licence to use UNCALLED for ReadUntil. This will grant you access to a private GitHub repository where you can install the ReadUntil API into your MinKNOW environment. You must install UNCALLED in the same way, and run the instance installed into MinKNOW environment (the above example assumes your MinKNOW python environment is located at `/opt/ONT/MinKNOW/ont-python`). We recommend that you try mapping fast5s with the MinKNOW installation via `uncalled map` before real-time enrichment, as runtime issues involving hdf5 libraries could come up if UNCALELD is not installed properly.

The command should be run once sequencing has begun and the "Channels Panel" has appeared. Running earlier may cause issues, including crashing MinKNOW and requiring you to restart your computer. Reads will not be ejected until the first mux scan has finished, which provides a ~1.5 minute window to start UNCALLED without missing any reads.

Arguments:

- `-x/--bwa-prefix` the prefix of the index to align to. Should be a BWA index that `uncalled index` was run on
- `-t/--threads` number of threads to use for mapping (default: 1)
- `-c/--max-chunks-proc` number of chunks to attempt mapping before giving up on a read (default: 10). One chunk is approximately one second of sequencing (~450bp).
- `--port` MinION device port. Use `uncalled list-ports` command to see all devices that have been plugged in since MinKNOW started.
- `--enrich` will *keep* reads that map to the reference if included
- `--deplete` will *eject* reads that map to the reference if included
- `--even` will only eject reads from even channels if included
- `--odd` will only eject reads from odd channels if included
- `--duration` expected duration of sequencing run in hours (default: 48)
- `--post-script` optional path to executable to run after analysis has finished. Useful to automatically basecall after sequencing is done, for example.
- `--post-time` buffer time to wait after the last read before running the post-script, in seconds (default: 60). Useful because MinKNOW's sequencing time is not always exact.

Note exactly one of `--deplete` or `--enrich` must be specified

## Output Format

Both `uncalled map` and `uncalled realtime` output to stdout in a format similar to [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md). Unmapped reads are output with reference-location-dependent fields replaced with \*s. Lines that begin with "#" are comments that useful for debugging.

Query coordinates, residue matches, and block lengths are estimated assuming 450bp sequenced per second. This estimate can be significantly off depending on the sequencing run. UNCALLED attempts to map a read as early as possible, so the "query end" field corresponds to the leftmost position where UNCALLED was able to confidently map the read. This differs from aligners such as [minimap2](https://github.com/lh3/minimap2), which attempt to map the full length of the read.

For real-time mapping, read lengths are estimated by how much signal UNCALLED recieved, which does not necessarily correspond to how much signal was actually sequenced.

Both modes include the following custom attributes in each PAF entry:

- `mt`: **map time**. Time in milliseconds it took to map the read.
- `ch`: **channel**. MinION channel that the read came from.
- `st`: **start sample**. Global _sequenicng_ start time of the read (in signal samples, 4000 samples/sec).

`uncalled realtime` also includes the following attributes:

- `ej`: **ejected**. Time that the eject signal was sent, in milliseconds since last chunk was received.
- `kp`: **kept**. Time that UNCALLED decided to keep the read, in milliseconds since last chunk was received.
- `en`: **ended**. Time that UNCALLED determined the read ended, in milliseconds since last chunk was received.
- `mx`: **mux scan**. Time that the read _would have_ been ejected, had it not have occured within a mux scan.
- `wt`: **wait time**. Time in milliseconds that the read was queued but was not actively being mapped, either due to thread delays or waiting for new chunks.

## Practical Considerations

For ReadUntil sequencing, the first decision to make is whether to perform **enrichment** or **depletion** (`--enrich` or `--deplete`). 
In enrichment mode, UNCALLED will eject a read if it *does not* map to the reference, meaning your target should be the reference. 
In depletion mode, UNCALLED will eject a read if it *does* map to the reference, meaning your target should be everything except your reference.

Note that enrichment necessitates a quick decision as to whether or not a read maps, since you want to eject a read as fast as possible. Usually ~95% of reads can be mapped within three seconds for highly non-reptetive references, so setting `-c/--max-chunks-proc` to `3` generally works well for enrichment. The default value of `10` works well for depletion.

UNCALLED currently does not support large (> ~100Mb) or highly repetitive references. 
The speed and mapping rate both progressively drop as references become larger and more repetitive. 
Bacterial genomes or small collections of bacterial genomes typically work well. 
Small segments of eukaryotic genomes can also be used, however the presence of any repetitvie elements will harm the performance. 
We hope to provide tools and/or guidelines for masking such references in the near future, and increasing the supported reference size and repeat tolerance is a long-term goal.

ReadUntil works best with longer reads. Maximize your read lengths for best results.

UNCALLED currently only supports reads sequenced with r9.4 chemistry.

## Undocumented Features

UNCALLED is a work in progress. Many parameters exist that are not documented here, but can be seen on the command line help information. Most users should leave these unchanged. They may be removed in future versions, or be replaced with hyperparameters that adjust the accuracy and speed of UNCALLED.

You may notice a `simulate` subcommand which I have also not documented. This feature is incomplete, may not be useful, and will likely be removed in future versions.

