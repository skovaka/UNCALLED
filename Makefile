LIBS=-lz -lm -lstdc++ -ldl
HDF5_LIB=-L/cm/shared/apps/hdf5/1.8.17/lib /cm/shared/apps/hdf5/1.8.17/lib/libhdf5.a
BOOST_INCLUDE=-I/cm/local/apps/boost/1.58.0/include
HDF5_INCLUDE=-I/cm/shared/apps/hdf5/1.8.17/include
CC=g++
CFLAGS=-Wall -std=c++11 -g
INCLUDE=-I./src/fast5/src -I./src ${BOOST_INCLUDE}



all: sigalign signal_extract

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(INCLUDE) $(HDF5_INCLUDE) 

sigalign: sigalign.o nano_fmi.o kmer_model.o
	$(CC) $(CFLAGS) nano_fmi.o kmer_model.o sigalign.o -o sigalign $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 

#signal_extract: signal_extract.o model_tools.o
#	$(CC) $(CFLAGS) model_tools.o signal_extract.o -o signal_extract $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 
	
# nano_bwt: nano_bwt.o
#
# 	$(CC) $(CFLAGS) nano_bwt.o -o nano_bwt $(LIBS) $(INCLUDE)

# nano_bwt.o: nano_bwt.cpp
# 	$(CC) $(CFLAGS) -c nano_bwt.cpp $(compile_opts)

clean:
	rm *.o
