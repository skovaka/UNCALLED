##Inputs
Reference genome and BED file containing a single insertion site

Note that the name field is required, and must be unique if you include multiple sites.
```
$ ls REF_FILES/
L.monocytogenes.fasta

$ grep ">" REF_FILES/L.monocytogenes.fasta
>Listeria_monocytogenes_complete_genome 2.992Mb

$ cat tgt.bed # Note BED start/end are the same to indicate insert coord
Listeria_monocytogenes_complete_genome  300000  300000  ins_test
```

##Create the reference and config file
Using a flank length of 100Kbp, which is much larger than you would ever want but it's easier for 

```
$ ~/code/UNCALLED/ins_scripts/create_ins_ref.sh \
    tgt.bed \                         # Insertion site
    100000 \                          # Flank length
    REF_FILES/L.monocytogenes.fasta \ # Reference genome
    REF_FILES/ins_test                # Output prefix

$ ls REF_FILES/
ins_test.fa  ins_test_flank.bed  ins_test.toml  L.monocytogenes.fasta  L.monocytogenes.fasta.fai

$ grep ">" REF_FILES/ins_test.fa # Use this as reference index
>ins_test::Listeria_monocytogenes_complete_genome:199999-299999
>ins_test::Listeria_monocytogenes_complete_genome:300001-400001

$ cat REF_FILES/ins_test.toml # Use this as input to "--conf" flag
[realtime]
fwd_tgts = [
        "ins_test::Listeria_monocytogenes_complete_genome:199999-299999",
]
rev_tgts = [
        "ins_test::Listeria_monocytogenes_complete_genome:300001-400001",
]
```

##Index the reference and run simulator

Note: you should probably do repeat masking before indexing for human

```
$ uncalled index REF_FILES/ins_test.fa
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.05 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.05 sec
Initializing parameter search
Computing default parameters
Writing default parameters
Done

$ uncalled sim REF_FILES/ins_test.fa $ctrl_fast5_pass --ctl-seqsum $ctrl_summary \
    --unc-seqsum $unc_summary --unc-paf $unc_pafs -t 1 -c 3 --sim-speed 0.25 \
    --selective --conf REF_FILES/ins_test.toml \
    > uncalled_out.paf 2> uncalled_err.txt
...

# Check reference name and strand for kept reads
$ grep kp: uncalled_out.paf | awk '{print $6"\t"$5}' | sort | uniq
ins_test::Listeria_monocytogenes_complete_genome:199999-299999  +
ins_test::Listeria_monocytogenes_complete_genome:300001-400001  -

# Check for ejected reads
$ grep ej: uncalled_out.paf | awk '{print $6"\t"$5}' | sort | uniq
*       *
ins_test::Listeria_monocytogenes_complete_genome:199999-299999  -
ins_test::Listeria_monocytogenes_complete_genome:300001-400001  +
```

The final commands show that all kept reads are upstream on the forward strand or downstream on the reverse, and all ejected reads are the oposite or unmapped
