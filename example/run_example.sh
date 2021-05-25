#!/bin/bash

#Expected output: f41a60f7-de4a-4b17-9f54-387e52d60b65	3086	186	251	+	Escherichia_coli_chromosome:2400000-2410000	10000	6826	6891	63	66	255

cd `dirname "${BASH_SOURCE[0]}"`

if [ -d index ]; then
    rm -r index
fi

mkdir index
uncalled index -o index/example_ref example_ref.fa
uncalled map -c 3 index/example_ref fast5_filename.txt
