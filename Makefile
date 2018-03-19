LIBS=-lz -lm -lstdc++ -ldl 
HDF5_LIB=-L/cm/shared/apps/hdf5/1.8.17/lib /cm/shared/apps/hdf5/1.8.17/lib/libhdf5.a
SDSL_LIB=-L/home-4/skovaka1@jhu.edu/software/sdsl-lite/lib  -lsdsl -ldivsufsort -ldivsufsort64
SDSL_INC=-I/home-4/skovaka1@jhu.edu/software/sdsl-lite/include
#BOOST_INCLUDE=-I/cm/local/apps/boost/1.58.0/include
HDF5_INCLUDE=-I/cm/shared/apps/hdf5/1.8.17/include
CC=g++
CFLAGS=-Wall -std=c++11 -O3
INCLUDE=-I./src/fast5/src -I./src/scrappie -I./src #${BOOST_INCLUDE}

all: uncalled dtw_test save_fmi test_fmi align_stats

#seed_tracker_test 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ $(SDSL_INC) $(INCLUDE) $(HDF5_INCLUDE) 

#sigalign: sigalign.o kmer_fmi.o kmer_model.o seed_graph.o
#	$(CC) $(CFLAGS) kmer_fmi.o kmer_model.o seed_graph.o sigalign.o -o sigalign $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS) 

uncalled: uncalled.o sdsl_fmi.o kmer_model.o seed_graph.o seed_tracker.o arg_parse.o basepairs.o range.o
	$(CC) $(CFLAGS) sdsl_fmi.o kmer_model.o seed_graph.o seed_tracker.o arg_parse.o uncalled.o basepairs.o range.o -o uncalled $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(SDSL_LIB) $(LIBS)

seed_tracker_test: seed_tracker_test.o kmer_model.o seed_tracker.o seed_graph.o kmer_fmi.o basepairs.o  range.o
	$(CC) $(CFLAGS) seed_tracker_test.o kmer_model.o kmer_fmi.o seed_tracker.o seed_graph.o basepairs.o range.o -o seed_tracker_test $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(SDSL_LIB)

dtw_test: dtw.o kmer_model.o arg_parse.o basepairs.o
	$(CC) $(CFLAGS) dtw.o kmer_model.o arg_parse.o basepairs.o -o dtw_test $(INCLUDE) $(HDF5_INCLUDE) $(HDF5_LIB) $(LIBS)

test_fmi: base_fmi.o sdsl_fmi.o test_fmi.o basepairs.o range.o
	$(CC) $(CFLAGS) base_fmi.o sdsl_fmi.o test_fmi.o basepairs.o range.o -o test_fmi  $(SDSL_LIB)

align_stats: base_fmi.o sdsl_fmi.o align_stats.o basepairs.o range.o
	$(CC) $(CFLAGS) base_fmi.o sdsl_fmi.o align_stats.o basepairs.o range.o -o align_stats  $(SDSL_LIB)

save_fmi: sdsl_fmi.o base_fmi.o save_fmi.o basepairs.o range.o
	$(CC) $(CFLAGS) sdsl_fmi.o base_fmi.o save_fmi.o basepairs.o range.o -o save_fmi  $(SDSL_LIB) $(LIBS)

#sdsl_fmi: sdsl_fmi.o basepairs.o range.o
#	$(CC) $(CFLAGS) sdsl_fmi.o basepairs.o range.o -o sdsl_fmi $(SDSL_LIB) $(LIBS)

#arg_parse_test: 
#	$(CC) $(CFLAGS) arg_parse.o -o arg_parse_test $(INCLUDE) $(LIBS)

clean:
	rm *.o
