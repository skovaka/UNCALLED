# UNCALLED: a Utility for Nanopore Current ALgnment to Large Expanses of DNA

UNCALLED is a rapid nanopore signal mapper intended for targeted sequencing via adaptive sampling/ReadUntil

Uncalled4 is an accurate end-to-end nanopore signal aligner for nucleotide modification detection and general-purpose signal analysis

Uncalled4 was initially developed as an extension of UNCALLED, but they are now two distinct projects. 

## Uncalled4: visualization and analysis of nanopore RNA and DNA signal alignments

Uncalled4 is a toolkit for alignment, analysis, and visualization of raw nanopore RNA and DNA signal. Uncalled4 features interactive alignment visualizations, an alignment algorithm guided by basecaller metadata, methods for training pore k-mer models, and epigenetic modification detection statistics. We also introduce a BAM signal alignment encoding which can be converted to and from [Nanopolish eventalign](https://github.com/jts/nanopolish)  and other signal alignment formats. For more information see https://github.com/skovaka/uncalled4

[![Uncalled4 logo](assets/img/logo4.png "Uncalled4 logo")](https://github.com/skovaka/uncalled4)

### UNCALLED: real-time targeted sequencing via raw signal mapping

UNCALLED enables targeted sequencing by directly aligning raw nanopore signal to a reference and selectively ejecting reads via [adaptive sampling](https://nanoporetech.com/resource-centre/adaptive-sampling-oxford-nanopore) (AKA ReadUntil).

For more information see the https://github.com/skovaka/UNCALLED, or [read the paper](https://www.nature.com/articles/s41587-020-0731-9)

