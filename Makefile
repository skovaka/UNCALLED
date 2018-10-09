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

UNCALLED_OBJS= kmer_model.o aligner.o seed_tracker.o arg_parse.o basepairs.o range.o

all: uncalled dtw_test save_fmi test_fmi align_stats detect_events

#uncalled_graph uncalled_forest 

#uncalled_graph.o: uncalled.cpp
#	$(CC) $(CFLAGS) -D ALN_TYPE=GRAPH_ALN -c -o $@ $^ $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 
#
#uncalled_forest.o: uncalled.cpp
#	$(CC) $(CFLAGS) -D ALN_TYPE=FOREST_ALN -c -o $@ $^ $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 
#

dtw_test.o: dtw_test.cpp dtw.hpp
	$(CC) $(CFLAGS) -c -o dtw_test.o dtw_test.cpp $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 

#uncalled_graph: $(UNCALLED_OBJS) uncalled_graph.o sdsl_fmi.o base_fmi.o graph_aligner.o
#	$(CC) $(CFLAGS) -D ALN_TYPE=GRAPH_ALN $(UNCALLED_OBJS) uncalled_graph.o sdsl_fmi.o base_fmi.o graph_aligner.o -o uncalled_graph $(HDF5_LIB) $(SDSL_LIB) $(LIBS)
#
#uncalled_forest: $(UNCALLED_OBJS) uncalled_forest.o sdsl_fmi.o base_fmi.o forest_aligner.o
#	$(CC) $(CFLAGS) -D ALN_TYPE=FOREST_ALN $(UNCALLED_OBJS) uncalled_forest.o sdsl_fmi.o base_fmi.o forest_aligner.o -o uncalled_forest $(HDF5_LIB) $(SDSL_LIB) $(LIBS)

detect_events: detect_events.o event_detector.o
	$(CC) $(CFLAGS) detect_events.o event_detector.o -o detect_events $(HDF5_LIB) $(LIBS)

uncalled: $(UNCALLED_OBJS) uncalled.o bwa_fmi.o leaf_arr_aligner.o event_detector.o
	$(CC) $(CFLAGS) $(UNCALLED_OBJS) uncalled.o bwa_fmi.o leaf_arr_aligner.o event_detector.o -o uncalled $(HDF5_LIB) $(BWA_LIB) $(LIBS)

dtw_test: dtw_test.o kmer_model.o arg_parse.o basepairs.o dtw.hpp event_detector.o
	$(CC) $(CFLAGS) dtw_test.o kmer_model.o arg_parse.o basepairs.o event_detector.o -o dtw_test $(INCLUDE) $(HDF5_LIB) $(LIBS)

test_fmi: base_fmi.o sdsl_fmi.o bwa_fmi.o test_fmi.o basepairs.o range.o
	$(CC) $(CFLAGS) base_fmi.o sdsl_fmi.o bwa_fmi.o test_fmi.o basepairs.o range.o -o test_fmi $(LIBS) $(SDSL_LIB) $(BWA_LIB)

align_stats: base_fmi.o sdsl_fmi.o bwa_fmi.o align_stats.o basepairs.o range.o
	$(CC) $(CFLAGS) base_fmi.o sdsl_fmi.o bwa_fmi.o align_stats.o basepairs.o range.o -o align_stats $(LIBS) $(SDSL_LIB) $(BWA_LIB)

save_fmi: sdsl_fmi.o base_fmi.o bwa_fmi.o save_fmi.o basepairs.o range.o
	$(CC) $(CFLAGS) sdsl_fmi.o base_fmi.o bwa_fmi.o save_fmi.o basepairs.o range.o -o save_fmi  $(SDSL_LIB) $(LIBS) $(BWA_LIB)

fast5_dump: fast5_dump.o 
	$(CC) $(CFLAGS) fast5_dump.o -o fast5_dump $(LIBS) $(INCLUDE) $(HDF5_LIB) $(HDF5_INCLUDE)

clean:
	rm *.o
