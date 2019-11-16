# Usage

```
> ./run_example.sh
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.00 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.00 sec
[main] Version: 0.7.17-r1194-dirty
[main] CMD: ../bwa/bwa index -p index/example_ref example_ref.fa
[main] Real time: 0.015 sec; CPU: 0.010 sec
Writing index/example_ref.uncl
Reading fast5 paths
Aligning
f41a60f7-de4a-4b17-9f54-387e52d60b65    3562    86      118     -       Escherichia_coli_chromosome:2400000-2410000    10000   6937    6975    38      39      255     ch:i:486
       st:i:257117    mt:f:22.596806
```

# Data Source

Example data obtained from https://github.com/LomanLab/mockcommunity

Read extracted from https://nanopore.s3.climb.ac.uk/Zymo-GridION-EVEN-BB-SN_signal.tar

Reference extracted from https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
