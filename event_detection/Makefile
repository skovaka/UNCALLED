LIBS=-lz -lm -lstdc++ -ldl
HDF5_LIB=-L/cm/shared/apps/hdf5/1.8.17/lib /usr/lib/libhdf5.a
HDF5_INCLUDE=-I/cm/shared/apps/hdf5/1.8.17/include
CC=gcc
CPP_C=g++
CPP_FLAGS=-Wall -std=c++11 -O3
C_FLAGS=-Wall -O3

all: evdt_test evdt_test_c

evdt_test.o: evdt_test.cpp
	$(CPP_C) $(CPP_FLAGS) -c -o $@ $^ $(HDF5_INCLUDE) 

event_detector.o: event_detector.cpp
	$(CPP_C) $(CPP_FLAGS) -c -o $@ $^

event_detection.o: event_detection.c
	$(CC) $(C_FLAGS) -c -o $@ $^ $(HDF5_INCLUDE) 

fast5_interface.o: fast5_interface.c
	$(CC) $(C_FLAGS) -c -o $@ $^ $(HDF5_INCLUDE) 

evdt_test_c: event_detection.o fast5_interface.o 
	$(CC) $(C_FLAGS) event_detection.o fast5_interface.o -o evdt_test_c $(HDF5_INCLUDE)  $(HDF5_LIB) $(LIBS) 

evdt_test: evdt_test.o fast5_interface.o event_detector.o
	$(CPP_C) $(CPP_FLAGS) evdt_test.o event_detector.o fast5_interface.o -o event_detect_test $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 

clean:
	rm *.o
