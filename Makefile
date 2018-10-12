LIBS=-lz -lm -lstdc++ -ldl 
HDF5_LIB=-L/home-4/skovaka1@jhu.edu/software/lib /home-4/skovaka1@jhu.edu/software/lib/libhdf5.a
SDSL_LIB=-L/home-4/skovaka1@jhu.edu/software/sdsl-lite/lib  -lsdsl -ldivsufsort -ldivsufsort64
BWA_LIB=-L/home-4/skovaka1@jhu.edu/code/nanopore_aligner/bwa /home-4/skovaka1@jhu.edu/code/nanopore_aligner/bwa/libbwa.a
SDSL_INC=-I/home-4/skovaka1@jhu.edu/software/sdsl-lite/include
#BOOST_INCLUDE=-I/cm/local/apps/boost/1.58.0/include
#HDF5_INCLUDE=-I/cm/shared/apps/hdf5/1.8.17/include
HDF5_INCLUDE=-I/home-4/skovaka1@jhu.edu/software/include
CC=g++
CFLAGS=-Wall -std=c++11 -O3
INCLUDE=-I./fast5/src #${BOOST_INCLUDE}

all: uncalled dtw_test align_stats detect_events

dtw_test.o: dtw_test.cpp dtw.hpp
	$(CC) $(CFLAGS) -c -o dtw_test.o dtw_test.cpp $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 

detect_events: detect_events.o event_detector.o
	$(CC) $(CFLAGS) detect_events.o event_detector.o -o detect_events $(HDF5_LIB) $(LIBS)

uncalled: uncalled.o kmer_model.o aligner.o seed_tracker.o arg_parse.o util.o range.o bwa_fmi.o event_detector.o
	$(CC) $(CFLAGS) uncalled.o kmer_model.o aligner.o seed_tracker.o arg_parse.o util.o range.o bwa_fmi.o event_detector.o -o uncalled $(HDF5_LIB) $(BWA_LIB) $(LIBS)

dtw_test: dtw_test.o kmer_model.o arg_parse.o util.o dtw.hpp event_detector.o
	$(CC) $(CFLAGS) dtw_test.o kmer_model.o arg_parse.o util.o event_detector.o -o dtw_test $(INCLUDE) $(HDF5_LIB) $(LIBS)

align_stats: bwa_fmi.o align_stats.o util.o range.o
	$(CC) $(CFLAGS) bwa_fmi.o align_stats.o util.o range.o -o align_stats $(LIBS) $(SDSL_LIB) $(BWA_LIB)

fast5_dump: fast5_dump.o 
	$(CC) $(CFLAGS) fast5_dump.o -o fast5_dump $(LIBS) $(INCLUDE) $(HDF5_LIB) $(HDF5_INCLUDE)

clean:
	rm *.o
