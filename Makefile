LIBS=-lstdc++ -lz -ldl 
HDF5_LIB=-L/home-4/skovaka1@jhu.edu/anaconda3/lib /home-4/skovaka1@jhu.edu/anaconda3/lib/libhdf5.a

BWA_LIB=-L/home-4/skovaka1@jhu.edu/code/nanopore_aligner/bwa /home-4/skovaka1@jhu.edu/code/nanopore_aligner/bwa/libbwa.a
HDF5_INCLUDE=-I/home-4/skovaka1@jhu.edu/software/include
CC=g++
CFLAGS=-Wall -std=c++11 -O3
INCLUDE=-I./fast5/src -I./pybind11/include #${BOOST_INCLUDE}

all: uncalled dtw_test self_align_ref detect_events

dtw_test.o: dtw_test.cpp dtw.hpp
	$(CC) $(CFLAGS) -c -o dtw_test.o dtw_test.cpp  $(INCLUDE) $(HDF5_INCLUDE) 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(INCLUDE) $(HDF5_INCLUDE) 

detect_events: detect_events.o event_detector.o
	$(CC) $(CFLAGS) detect_events.o event_detector.o -o detect_events $(HDF5_LIB) $(LIBS)

uncalled: uncalled.o kmer_model.o aligner.o seed_tracker.o arg_parse.o range.o bwa_fmi.o event_detector.o
	$(CC) $(CFLAGS) uncalled.o kmer_model.o aligner.o seed_tracker.o arg_parse.o range.o bwa_fmi.o event_detector.o -o uncalled $(HDF5_LIB) $(BWA_LIB) $(LIBS)

dtw_test: dtw_test.o kmer_model.o arg_parse.o dtw.hpp event_detector.o
	$(CC) $(CFLAGS) dtw_test.o kmer_model.o arg_parse.o event_detector.o -o dtw_test $(INCLUDE) $(HDF5_LIB) $(LIBS)

self_align_ref: bwa_fmi.o self_align_ref.o range.o arg_parse.o
	$(CC) $(CFLAGS) bwa_fmi.o self_align_ref.o range.o arg_parse.o -o self_align_ref $(LIBS) $(BWA_LIB)

fast5_dump: fast5_dump.o 
	$(CC) $(CFLAGS) fast5_dump.o -o fast5_dump $(LIBS) $(INCLUDE) $(HDF5_LIB) $(HDF5_INCLUDE)

clean:
	rm *.o
