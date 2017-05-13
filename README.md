# ont_sig_align
Nanopore signal aligner


# Prerequisites:

- HDF5:
```
cd src
wget https://support.hdfgroup.org/ftp/HDF5/current18/bin/linux-centos7-x86_64-gcc485/hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared.tar.gz
tar -xzvf hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared.tar.gz
```

- Boost:
```
# make sure boost is installed
sudo apt-get install libboost1.58-all-dev
# find where the boost libraries are located on your computer, set to $BOOSTDIR
cd src
bcp --boost=$BOOSTDIR math/distributions boost/
```
