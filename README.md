# ont_sig_align
Nanopore signal aligner



# Installation

## Prerequisites:

-Fast5:
```
cd src
git clone https://github.com/mateidavid/fast5
make
```

- HDF5:
Should be included on the MARCC server at location specified in Makefile. Otherwise:
```
cd src
wget https://support.hdfgroup.org/ftp/HDF5/current18/bin/linux-centos7-x86_64-gcc485/hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared.tar.gz
tar -xzvf hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared.tar.gz
```

- Boost:
Should be included on the MARCC server at location specified in Makefile. Otherwise:
```
# make sure boost is installed
sudo apt-get install libboost1.58-all-dev
# find where the boost libraries are located on your computer, set to $BOOSTDIR
cd src
bcp --boost=$BOOSTDIR math/distributions boost/
```

## Compilation
make sure `BOOST_INCLUDE`, `HDF5_INCLUDE` and `HDF5_LIB` are set accordingly in Makefile

```
make
```


# running

```
./sigalign REFERENCE.FASTA MODEL TALLY_DISTANCE READ1.FAST5 READ2.FAST5 ...
```

expected output:
```
READ1.FAST5: <x> matches on fwd strand, <y> matches on rev strand
READ2.FAST5: <x> matches on fwd strand, <y> matches on rev strand
```
