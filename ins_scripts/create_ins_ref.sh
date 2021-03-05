#!/bin/bash

INS_SITE_BED=$1
FLANK_LEN=$2
GENOME_FASTA=$3
OUT_PREFIX=$4

FAIDX=${GENOME_FASTA}.fai
FLANK_BED=${OUT_PREFIX}_flank.bed

OUT_FASTA=${OUT_PREFIX}.fa
OUT_TOML=${OUT_PREFIX}.toml

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#Create reference index if not already present
if [ ! -f ${FAIDX} ]; then
    samtools faidx ${GENOME_FASTA}
fi

#Create flanking sequnce bed
bedtools flank -i ${INS_SITE_BED} -b ${FLANK_LEN} -g ${FAIDX} > ${FLANK_BED}

#Pull flanking sequence from fasta
bedtools getfasta -name -bed ${FLANK_BED} -fi ${GENOME_FASTA} > ${OUT_FASTA}

#Write UNCALLED config file to specify upstream/downstream sequences
python3 ${SCRIPT_DIR}/write_conf.py ${OUT_FASTA} > ${OUT_TOML}
