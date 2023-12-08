# UNCALLED: a Utility for Nanopore Current ALgnment to Large Expanses of DNA

The UNCALLED repository contains two related tools for nanopore signal alignment/mapping: UNCALLED for signal mapping and adaptive sampling, and Uncalled4 for accurate signal alignment, analysis, and visualization.

## Uncalled4: visualization and analysis of nanopore RNA and DNA signal alignments

We are actively developing Uncalled4: a toolkit for alignment, analysis, and visualization of raw nanopore RNA and DNA signal. Uncalled4 features interactive alignment visualizations, an alignment algorithm guided by basecaller metadata, methods for training pore k-mer models, and epigenetic modification detection statistics. This was initially developed as an extension of UNCALLED (see below), but we now plan to migrate it to a dedicated repository soon. In the meantime it is currently availible in the [uncalled4 branch](https://github.com/skovaka/UNCALLED/tree/uncalled4) of the UNCALLED repository.

[![Uncalled4 logo](assets/img/logo4.png "Uncalled4 logo")](https://github.com/skovaka/UNCALLED/tree/uncalled4)

### UNCALLED: real-time targeted sequencing via raw signal mapping

UNCALLED enables targeted sequencing by directly aligning raw nanopore signal to a reference and selectively ejecting reads via [adaptive sampling](https://nanoporetech.com/resource-centre/adaptive-sampling-oxford-nanopore) (AKA ReadUntil). This is a distinct algorithm from Uncalled4 (above), focused on rapid seed mapping rather than accurate end-to-end alignment.

For more information see the [main branch](https://github.com/skovaka/UNCALLED), or [read the paper](https://www.nature.com/articles/s41587-020-0731-9)

