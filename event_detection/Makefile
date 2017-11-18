LIBS=-lz -lm -lstdc++ -ldl
HDF5_LIB=-L/cm/shared/apps/hdf5/1.8.17/lib /cm/shared/apps/hdf5/1.8.17/lib/libhdf5.a
HDF5_INCLUDE=-I/cm/shared/apps/hdf5/1.8.17/include
CC=gcc
CPP_C=g++
CPP_FLAGS=-Wall -std=c++11 -O3
C_FLAGS=-Wall -O3

all: event_detect_test

event_detect_test.o: event_detect_test.cpp
	$(CPP_C) $(CPP_FLAGS) -c -o $@ $^ $(HDF5_INCLUDE) 

evenit_detector.o: event_detector.cpp
	$(CPP_C) $(CPP_FLAGS) -c -o $@ $^

fast5_interface.o: fast5_interface.c
	$(CC) $(C_FLAGS) -c -o $@ $^ $(HDF5_INCLUDE) 

event_detect_test: event_detect_test.o fast5_interface.o event_detector.o
	$(CPP_C) $(CPP_FLAGS) event_detect_test.o event_detector.o fast5_interface.o -o event_detect_test $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 

#signal_extract: signal_extract.o model_tools.o
#	$(CC) $(CFLAGS) model_tools.o signal_extract.o -o signal_extract $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 
	
# nano_bwt: nano_bwt.o
#
# 	$(CC) $(CFLAGS) nano_bwt.o -o nano_bwt $(LIBS) $(INCLUDE)

# nano_bwt.o: nano_bwt.cpp
# 	$(CC) $(CFLAGS) -c nano_bwt.cpp $(compile_opts)

clean:
	rm *.o
