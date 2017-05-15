LIBS=-lsz -lz -lm -lstdc++ -ldl
HDF5_LIB=-L./src/hdf5/lib ./src/hdf5/lib/libhdf5.a
HDF5_INCLUDE=-I./src/hdf5/include
CC=g++
CFLAGS=-Wall -std=c++11 -O3
INCLUDE=-I./src/fast5/src -I./src -I./src/boost



all: sigalign

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(INCLUDE) $(HDF5_INCLUDE) 

sigalign: sigalign.o nano_bwt.o model_tools.o
	$(CC) $(CFLAGS) nano_bwt.o model_tools.o sigalign.o -o sigalign $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 

# nano_bwt: nano_bwt.o
# 	$(CC) $(CFLAGS) nano_bwt.o -o nano_bwt $(LIBS) $(INCLUDE)

# nano_bwt.o: nano_bwt.cpp
# 	$(CC) $(CFLAGS) -c nano_bwt.cpp $(compile_opts)

clean:
	rm *.o
