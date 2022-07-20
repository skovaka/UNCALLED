# Simulator Scripts

This directory contains two scripts which can be used to interpet the output of `uncalled sim`.

We understand these scripts may not be the most user friendly. We will work towards improving them and better intergrate them with the simulator in the future.

## `est_genome_yield.py`

**Example:**
```
> sim_scripts/est_genome_yield.py -u uncalled_out.paf --enrich -x E.coli -m mm2.paf -s sequencing_summary.txt --sim-speed 0.25

unc_on_bp       150.678033
unc_total_bp    6094.559395
cnt_on_bp       33.145022
cnt_total_bp    8271.651331
```

This is designed to be used in the context of enriching/depleting for whole genomes or chromosomes.

Arguments:
- `-u/--uncalled-fname`: Simulator output PAF file
- `-s/--seq_sum`: Control sequencing summary
- `-m/--minimap-fname`: Minimap2 PAF file of the control reads aligned to a reference containing the target (or off-target, in the case of depletion) sequences
- `-x/--bwa-prefix`: BWA reference used during the simulation
- `--deplete/--enrich`: same as option used in simulation
- `-t/--sim-speed`: Speed that the simulator was run at in the range (0.0, 1.0]

## `est_bed_yield.py`

**Example:**
```
> sim_scripts/est_bed_yield.py -u uncalled_out.paf -c ctl_coverage.bed -s sequencing_summary.txt -t 0.25 

unc_on_bp       150.678033
unc_total_bp    6094.559395
cnt_on_bp       33.145022
cnt_total_bp    8271.651331

```

This is designed to be used when the targets are subsequences of a larger reference, for example a set of genes. This requires `bedtools intersect -bed -a control_alns.bam -b targets.bed` to be run prior, where `control_alns.bam` is minimap2 alignments of the basecalled control reads to the full reference, and `targets.bed` are the targeted regions.

Arguments:
- `-u/--uncalled-fname`: Simulator output PAF file
- `-c/--cov-fname`: BED file of control read coverage. Should be output from 'bedools intersect' of control read alignments and the target region(s)
- `-s/--seq-sum`: Control sequencing summary
- -t/--sim-speed: Speed that the simulator was run at
