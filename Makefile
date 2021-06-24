CC=gcc
CFLAGS=-Wall -std=c++11 -g -fPIC -O3 $(FLAGS)

LIBHDF5=./submods/hdf5/lib/libhdf5.a
HDF5_LIB=-L./submods/hdf5/lib $(LIBHDF5)
HDF5_INCLUDE=-I./submods/hdf5/include

LIBBWA=./submods/bwa/libbwa.a
BWA_LIB=-L./submods/bwa $(LIBBWA)
BWA_INCLUDE=-I./submods/bwa

LIBS=$(HDF5_LIB) $(BWA_LIB) -lstdc++ -lz -ldl -pthread -lm 
INCLUDE=-I submods/ -I submods/toml11 -I submods/pybind11/include -I submods/pdqsort $(HDF5_INCLUDE) $(BWA_INCLUDE)

SRC=src
BUILD=build
BIN=bin
LIB=lib
#INCLUDE=include

_COMMON_OBJS=mapper.o seed_tracker.o range.o event_detector.o normalizer.o read_buffer.o fast5_reader.o event_profiler.o paf.o pore_model.o

_MAP_ORD_OBJS=$(_COMMON_OBJS) realtime_pool.o map_pool_ord.o uncalled_map_ord.o 
_MAP_OBJS=$(_COMMON_OBJS) map_pool.o uncalled_map.o 
_SIM_OBJS=$(_COMMON_OBJS) realtime_pool.o simulator.o uncalled_sim.o 
_DTW_OBJS=dtw_test.o fast5_reader.o read_buffer.o normalizer.o event_detector.o range.o event_profiler.o pore_model.o

_ALL_OBJS=$(_COMMON_OBJS) realtime_pool.o map_pool.o uncalled_map.o uncalled_map_ord.o simulator.o uncalled_sim.o dtw_test.o

MAP_OBJS = $(patsubst %, $(BUILD)/%, $(_MAP_OBJS))
MAP_ORD_OBJS = $(patsubst %, $(BUILD)/%, $(_MAP_ORD_OBJS))
SIM_OBJS = $(patsubst %, $(BUILD)/%, $(_SIM_OBJS))
DTW_OBJS = $(patsubst %, $(BUILD)/%, $(_DTW_OBJS))
ALL_OBJS = $(patsubst %, $(BUILD)/%, $(_ALL_OBJS))

DEPENDS := $(patsubst %.o, %.d, $(ALL_OBJS))

MAP_BIN = $(BIN)/uncalled_map
MAP_ORD_BIN = $(BIN)/uncalled_map_ord
SIM_BIN = $(BIN)/uncalled_sim
DTW_BIN = $(BIN)/dtw_test

all: dirs $(MAP_BIN) $(MAP_ORD_BIN) $(SIM_BIN) $(DTW_BIN)

#$(BIN)/%.o:src/%.c
#	$(CC) -c $< -o $@

$(MAP_BIN): $(MAP_OBJS) $(LIBHDF5) $(LIBBWA)
	$(CC) $(CFLAGS) $(MAP_OBJS) -o $@ $(LIBS)

$(MAP_ORD_BIN): $(MAP_ORD_OBJS) $(LIBHDF5) $(LIBBWA)
	$(CC) $(CFLAGS) $(MAP_ORD_OBJS) -o $@ $(LIBS)

$(SIM_BIN): $(SIM_OBJS) $(LIBHDF5) $(LIBBWA)
	$(CC) $(CFLAGS) $(SIM_OBJS) -o $@ $(LIBS)

$(DTW_BIN): $(DTW_OBJS) $(LIBHDF5) $(LIBBWA)
	$(CC) $(CFLAGS) $(DTW_OBJS) -o $@ $(LIBS)
	
#inspired by https://github.com/jts/nanopolish/blob/master/Makefile
$(LIBHDF5):
	cd submods/hdf5 && \
		./configure --enable-threadsafe --disable-hl --prefix=`pwd` --enable-shared=no --with-pic=yes || exit 255
	make -j 8 -C submods/hdf5 && make -C submods/hdf5 install
#--CFLAGS="-fPIC"

$(LIBBWA):
	make -C submods/bwa -f ../../src/Makefile_bwa

-include $(DEPENDS)

$(BUILD)/%.o: $(SRC)/%.cpp $(LIBHDF5) $(LIBBWA)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@ $(INCLUDE)

DIRS: $(BIN) $(BUILD)

.PHONY: dirs
dirs: $(BIN)/ $(BUILD)/

$(BIN)/:
	mkdir -p $@

$(BUILD)/:
	mkdir -p $@

clean:
	rm -rf $(BUILD)
