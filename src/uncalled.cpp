#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "self_align_ref.hpp"
#include "mapper.hpp"
#include "fast5_pool.hpp"
#include "channel_pool.hpp"
#include "chunk_pool.hpp"
#include "fast5_reader.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(align, m) {
    m.doc() = "UNCALLED";

    py::class_<MapperParams>(m, "MapperParams")
        .def(py::init<const std::string&, //index_fname
                      const std::string&, //model_fname
                      u32,   //seed_len
                      u32,   //min_aln_len
                      u32,   //min_rep_len
                      u32,   //max_rep_copy
                      u32,   //max_consec_stay
                      u32,   //max_paths
                      u32,   //max_events_proc
                      u32,   //max_chunks_proc
                      u32,   //evt_buffer_len
                      u32,   //evt_winlen1
                      u32,   //evt_winlen2
                      u16,   //evt_batch_size
                      float, //evt_timeout
                      float, //evt_thresh1
                      float, //evt_thresh2
                      float, //evt_peak_height
                      float, //evt_min_mean
                      float, //evt_max_mean
                      float, //max_stay_frac
                      float, //min_seed_prob
                      float, //min_mean_conf
                      float  //min_top_conf
                      >());

    py::class_<Mapper>(m, "Mapper")
        .def(py::init<MapperParams &, u16>())
        .def("map_fast5", &Mapper::map_fast5);

    py::class_<ReadLoc>(m, "ReadLoc")
        .def(py::init())
        .def("str", &ReadLoc::str)
        .def("is_valid", &ReadLoc::is_valid)
        .def("get_channel", &ReadLoc::get_channel)
        .def("get_number", &ReadLoc::get_number)
        .def("set_unblocked", &ReadLoc::set_unblocked)
        .def("get_ref", &ReadLoc::get_ref);

    py::class_<Fast5Pool>(m, "Fast5Pool")
        .def(py::init<MapperParams &, u16, u32>())
        .def("add_fast5s", &Fast5Pool::add_fast5s)
        .def("update", &Fast5Pool::update)
        .def("all_finished", &Fast5Pool::all_finished)
        .def("stop_all", &Fast5Pool::stop_all); 

    py::class_<ChannelPool>(m, "ChannelPool")
        .def(py::init<MapperParams &, u16, u16>()) 
        .def("add_fast5s", &ChannelPool::add_fast5s)
        .def("update", &ChannelPool::update)
        .def("all_finished", &ChannelPool::all_finished)
        .def("stop_all", &ChannelPool::stop_all); 

    py::class_<ChunkPool>(m, "ChunkPool")
        .def(py::init<MapperParams &, u16, u16>()) 
        .def("update", &ChunkPool::update)
        .def("all_finished", &ChunkPool::all_finished)
        .def("stop_all", &ChunkPool::stop_all)
        .def("add_chunk", &ChunkPool::add_chunk)
        .def("end_read", &ChunkPool::end_read);

    py::class_<Chunk>(m, "Chunk")
        .def(py::init<const std::string&, //_id, 
                      u32, //_number, 
                      u64, //_chunk_start_sample, 
                      const std::vector<float>, //&_raw_data, 
                      u32, //raw_st, 
                      u32  //raw_len
                     >())
        .def("get_number", &Chunk::get_number)
        .def("size", &Chunk::size);

    py::class_<ChunkSim>(m, "ChunkSim")
        .def(py::init<u32, //max_loaded,
                      u32, //num_chs
                      u16, //chunk_len
                      float //speed
                     >())
        .def("add_files", &ChunkSim::add_files)
        .def("add_reads", &ChunkSim::add_reads)
        .def("get_read_chunks", &ChunkSim::get_read_chunks)
        .def("stop_receiving_read", &ChunkSim::stop_receiving_read)
        .def("unblock", &ChunkSim::unblock)
        .def("set_time", &ChunkSim::set_time)
        .def("start", &ChunkSim::start)
        .def("is_running", &ChunkSim::is_running);
    
    m.def("self_align", &self_align);
}

