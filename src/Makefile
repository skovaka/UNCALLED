LIBS=-lstdc++ -lz -ldl -pthread

#HDF5_LIB=-L/home-4/skovaka1@jhu.edu/anaconda3/lib /home-4/skovaka1@jhu.edu/anaconda3/lib/libhdf5.a
#BWA_LIB=-L/home-4/skovaka1@jhu.edu/code/nanopore_aligner/bwa /home-4/skovaka1@jhu.edu/code/nanopore_aligner/bwa/libbwa.a

HDF5_LIB=-lhdf5
BWA_LIB=-L/home/skovaka/Dropbox/code/jhu/UNCALLED/bwa /home/skovaka/Dropbox/code/jhu/UNCALLED/bwa/libbwa.a
HDF5_INCLUDE=-I/home-4/skovaka1@jhu.edu/software/include
CC=g++
CFLAGS=-Wall -std=c++11 -O3
INCLUDE=-I../fast5/include -I../pybind11/include -I../ -I../pdqsort #${BOOST_INCLUDE}

all: simulator_test map_test
#uncalled dtw_test self_align_ref detect_events

#dtw_test.o: dtw_test.cpp dtw.hpp
#	$(CC) $(CFLAGS) -c -o dtw_test.o dtw_test.cpp  $(INCLUDE) $(HDF5_INCLUDE) 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(INCLUDE) $(HDF5_INCLUDE) 

map_test: map_test.o kmer_model.o mapper.o seed_tracker.o range.o bwa_fmi.o event_detector.o chunk_pool.o normalizer.o chunk.o params.o read_buffer.o fast5_pool.o
	$(CC) $(CFLAGS) map_test.o kmer_model.o mapper.o seed_tracker.o range.o bwa_fmi.o event_detector.o chunk_pool.o normalizer.o chunk.o params.o read_buffer.o fast5_pool.o -o map_test $(HDF5_LIB) $(BWA_LIB) $(LIBS)

simulator_test: simulator_test.o kmer_model.o mapper.o seed_tracker.o range.o bwa_fmi.o event_detector.o chunk_pool.o simulator.o normalizer.o chunk.o params.o read_buffer.o 
	$(CC) $(CFLAGS) simulator_test.o kmer_model.o mapper.o seed_tracker.o range.o bwa_fmi.o event_detector.o chunk_pool.o simulator.o normalizer.o chunk.o params.o read_buffer.o -o simulator_test $(HDF5_LIB) $(BWA_LIB) $(LIBS)

#uncalled: uncalled.o kmer_model.o mapper.o seed_tracker.o arg_parse.o range.o bwa_fmi.o event_detector.o chunk_pool.o fast5_reader.o
#	$(CC) $(CFLAGS) uncalled.o kmer_model.o mapper.o seed_tracker.o arg_parse.o range.o bwa_fmi.o event_detector.o chunk_pool.o fast5_reader.o -o uncalled $(HDF5_LIB) $(BWA_LIB) $(LIBS)
#
#dtw_test: dtw_test.o kmer_model.o arg_parse.o dtw.hpp event_detector.o
#	$(CC) $(CFLAGS) dtw_test.o kmer_model.o arg_parse.o event_detector.o -o dtw_test $(INCLUDE) $(HDF5_LIB) $(LIBS)
#
#self_align_ref: bwa_fmi.o self_align_ref.o range.o arg_parse.o
#	$(CC) $(CFLAGS) bwa_fmi.o self_align_ref.o range.o arg_parse.o -o self_align_ref $(LIBS) $(BWA_LIB)
#
#fast5_dump: fast5_dump.o 
#	$(CC) $(CFLAGS) fast5_dump.o -o fast5_dump $(LIBS) $(INCLUDE) $(HDF5_LIB) $(HDF5_INCLUDE)

clean:
	rm *.o
