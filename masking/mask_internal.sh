print_usage() {
   echo "Error: $1" 1>&2
    echo 1>&2
    echo "Usage: mask_internal.sh <reference> <k> <iters> <out_prefix>" 1>&2
    echo 1>&2
    echo "reference  - fasta file of the reference to mask" 1>&2
    echo "k          - k-mer length" 1>&2
    echo "iters      - number of masking iterations to run" 1>&2
    echo "out_prefix - output prefix for all output files" 1>&2
}

reference=$1
k=$2
iters=$3
prefix=$4

i=0

if [ ! -f "${reference}" ]; then
    print_usage "please provide reference fasta"
    exit 1
fi

if [ $k -lt 1 ]; then
    print_usage "k-mer length must be number greater than 0"
    exit 1
fi

if [ $iters -lt 1 ]; then
    print_usage "iterations must be number greater than 0"
    exit 1
fi

if [ ! -d `dirname ${prefix}` ]; then
    print_usage "directory \"`dirname ${prefix}`\" does not exist"
    exit 1
fi

if [ ! -f "${prefix}mask0.fa" ]; then
    ln ${reference} "${prefix}mask0.fa"
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


for j in `seq $iters`; do
    jellyfish count -o "${prefix}mask${i}.jf" -m $k -s 1000000 "${prefix}mask${i}.fa"
    c=`jellyfish stats "${prefix}mask${i}.jf" | grep Max_count | awk '{print $2}'`
    kmer=`jellyfish dump "${prefix}mask${i}.jf" -L ${c} | tail -n 1`
 
    printf "Iteration $i: " 1>&2
    ${SCRIPT_DIR}/mask_kmers.py ${prefix}mask${i}.fa -k $kmer > "${prefix}mask${j}.fa" 
    i=${j}
done
