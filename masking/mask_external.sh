print_usage() {
   echo "Error: $1" 1>&2
    echo 1>&2
    echo "Usage: mask_external.sh <full_bowtie> <target_fa> <min_len> <min_copy> <threads> <out_prefix>" 1>&2
    echo 1>&2
    echo "full_bowtie - bowtie prefix of index search for repeats of sequences in the target" 1>&2   
    echo "target_fa   - FASTA file of the target reference" 1>&2
    echo "min_len     - minimum repeat length" 1>&2
    echo "min_copy    - minimum repeat copy number" 1>&2
    echo "threads     - number of threads to use for bowtie" 1>&2 
    echo "out_prefix  - output prefix for all output files (default: working directory)" 1>&2
}

bt_prefix=$1
tgt_fa=$2
k=$3
min_copy=$4
threads=$5
out_prefix=$6

if [ ! -f "${bt_prefix}.1.ebwt" ]; then
    print_usage "${bt_prefix} does not appear to be a bowtie index prefix"
    exit 1
fi

if [ ! -f "${tgt_fa}" ]; then
    print_usage "target FASTA ${tgt_fa} does not exist"
    exit 1
fi

int_re='^[0-9]+$'

if [[ ! $k =~ $int_re ]] || [ $k -lt 2 ]; then
    print_usage "min_len must be integer greater than 1"
    exit 1
fi

if [[ ! $min_copy =~ $int_re ]] || [ $min_copy -lt 1 ]; then
    print_usage "min_copy must be integer greater than 0"
    exit 1
fi

if [[ ! $threads =~ $int_re ]] || [ $threads -lt 1 ]; then
    print_usage "threads must be integer greater than 0"
    exit 1
fi

if [ ! -d `dirname "${out_prefix}"` ]; then
    print_usage "directory \"`dirname ${out_prefix}`\" does not exist"
    exit 1
fi

ref_bed=${ref_prefix}.bed
compl=${out_prefix}ref_compl
windows=${out_prefix}windows
tgt_sizes=${out_prefix}tgt.sizes


samtools faidx ${tgt_fa}
awk '{print $1"\t"$2}' ${tgt_fa}.fai > ${tgt_sizes}
bedtools makewindows -g ${tgt_sizes} -w ${k} -s 1 | awk -v k=$k '$3-$2 == k' > ${windows}.bed
bedtools getfasta -fi ${tgt_fa} -bed ${windows}.bed > ${windows}.fa

rep_txt=${out_prefix}reps.txt
rep_bed=${out_prefix}reps.bed

bowtie -p ${threads} -fa -v 0 --suppress 2,3,4,5,6,7 ${bt_prefix} ${windows}.fa 2> /dev/null\
    | uniq -c \
    | awk '$1 > 1 {print $2"\t"$1}' \
    | sort -n -k 2,2 > ${rep_txt}

awk -F '[\t:-]' '{print $1":"$2"-"$3"\t"$4"\t"$5"\t"$6}' ${rep_txt} > ${rep_bed}

mask_bed=${out_prefix}reps_m${min_copy}.bed
mask_fa=${out_prefix}masked${min_copy}.fa

awk -v min_copy=$min_copy '$4 > min_copy' ${rep_bed} \
   | bedtools sort | bedtools merge > ${mask_bed}

bedtools maskfasta -fi ${tgt_fa} -bed ${mask_bed} -fo $mask_fa

masked_bp=`awk '{s+=($3-$2)} END {print s}' ${mask_bed}`

echo "Masked $masked_bp basepairs"
