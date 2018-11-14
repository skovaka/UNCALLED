LIBS=-lz -lm -lstdc++ -ldl -lhdf5 -lbwa
CC=g++
CFLAGS=-Wall -fPIC -std=c++11 -O3
INCLUDE=-I../fast5/include -I../pybind11/include -I../
LIB_DIRS=-L../bwa

all: normalize_test

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(INCLUDE)

normalize_test: normalizer.o event_detector.o kmer_model.o normalize_test.o
	$(CC) $(CFLAGS) normalizer.o event_detector.o kmer_model.o normalize_test.o -o normalize_test $(LIB_DIRS) $(LIBS)

clean:
	rm *.o
