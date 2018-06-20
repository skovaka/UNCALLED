#
# Part of: https://github.com/mateidavid/fast5
#
# (c) 2017: Matei David, Ontario Institute for Cancer Research
# MIT License
#

.SUFFIXES:
MAKEFLAGS += -r
SHELL := /bin/bash
.DELETE_ON_ERROR:
.PHONY: all help list clean check_hdf5

HDF5_DIR ?= /cm/shared/apps/hdf5/1.8.17/
HDF5_INCLUDE_DIR ?= ${HDF5_DIR}/include
HDF5_LIB_DIR ?= ${HDF5_DIR}/lib
HDF5_LIB ?= hdf5
TCLAP_DIR ?= tclap
HPPTOOLS_DIR ?= hpptools

TARGETS = f5ls f5ls-full hdf5-mod f5-mod
EXTRA_TARGETS = f5dump f5pack
HPP_FILES := fast5.hpp hdf5_tools.hpp Huffman_Packer.hpp Bit_Packer.hpp

CXXFLAGS := -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -Wpedantic
CPPFLAGS := -isystem ${HDF5_INCLUDE_DIR}
EXTRA_CPPFLAGS := -isystem ${TCLAP_DIR}/include -I ${HPPTOOLS_DIR}/include
LDFLAGS := -L${HDF5_LIB_DIR} -Wl,--rpath=${HDF5_LIB_DIR} -l${HDF5_LIB} -lpthread -lz -ldl

default: ${TARGETS}

all: default ${EXTRA_TARGETS}

print-%:
	@echo '$*=$($*)'

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

list:
	@echo "TARGETS=${TARGETS}"

clean:
	rm -rf ${TARGETS} ${EXTRA_TARGETS}

check_hdf5:
	@[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }
	@[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }

check_tclap:
	@[ -f "${TCLAP_DIR}/include/tclap/CmdLine.h" ] || { echo "TCLAP not found; get it from https://github.com/mateidavid/tclap.git" >&2; exit 1; }

check_hpptools:
	@[ -f "${HPPTOOLS_DIR}/include/alg.hpp" ] || { echo "HPPTOOLS not found; get it from https://github.com/mateidavid/hpptools.git" >&2; exit 1; }

%: %.cpp ${HPP_FILES} | check_hdf5
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -o $@ $< ${LDFLAGS}

f5dump: f5dump.cpp ${HPP_FILES} | check_hdf5 check_tclap check_hpptools
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${EXTRA_CPPFLAGS} -o $@ $< ${LDFLAGS}

f5pack: f5pack.cpp ${HPP_FILES} File_Packer.hpp | check_hdf5 check_tclap check_hpptools
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${EXTRA_CPPFLAGS} -o $@ $< ${LDFLAGS}
