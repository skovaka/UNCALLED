# Repeat Masking

UNCALLED is currently slower and less accurate on references with many repeats and low-complexity sequences. We have included two scripts designed to mask the types of repeats which cause problems for UNCALLED. 

We recommend using these scripts if any intergenic eukaryotic DNA is included in your target reference. They are not necessary for most bacterial DNA references.

## `mask_internal.sh`: high-frequency "internal" k-mer masking

This script iteratively masks high-frequency k-mers within the reference to be input to UNCALLED. Each iteration it uses jellyfish to identify the highest frequency k-mer and then masks it. An iterative approach is used because masking one k-mer will affect the frequency of overlapping k-mers, so the k-mer counts must be re-computed every time.

**Requirements:**
- python2 or 3
- jellyfish

```
Usage: mask_internal.sh <reference> <k> <iters> <out_prefix>

reference  - fasta file of the reference to mask
k          - k-mer length
iters      - number of masking iterations to run
out_prefix - output prefix for all output files
```

This script improves the true positive mapping rate and speed, at the expense of increasing the false positive rate. We recommend running this first, then using `mask_external.sh` to reduce the false positive rate.

## `mask_external.sh`: long exact repeat masking

This script masks repeats that may occur between the target sequence a reference containing all the DNA you expect to be sequencing, including off-target sequences. For example, if the target is a set of human genes, the full reference should be the human genome. It is based on exact bowtie alignments of segments of the target reference.

**Requirements:**
- python2 or 3
- bedtools
- samtools
- bowtie (full reference must be indexed with bowtie_build)

```
Usage: mask_external.sh <full_bowtie> <target_fa> <min_len> <min_copy> <threads> <out_prefix>

full_bowtie - bowtie prefix of index search for repeats of sequences in the target
target_fa   - FASTA file of the target reference
min_len     - minimum repeat length
min_copy    - minimum repeat copy number
threads     - number of threads to use for bowtie
out_prefix  - output prefix for all output files (default: working directory)
```

This script reduces the false positive mapping rate by preventing off-target DNA from falsely mapping to a repeat in the target reference. It also can increase the true positive rate and speed.

## Recommendations

We have found that 30 iterations if internal masking with k-mer length 10 works well, followed by expernal masking with `min_len` of 50bp and `min_copy` of 5. In other words, simply run:

```
mask_internal.sh <target_reference> 10 30 <out>
mask_external.sh <full_reference> <out>mask30.fa 50 5 <threads> <out>mask30_
```

## Future work

We plan to further refine these scripts and integrate them with UNCALLED's indexing procedure in the future.
