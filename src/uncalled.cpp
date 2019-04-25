#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "self_align_ref.hpp"
#include "mapper.hpp"
#include "simulator.hpp"
//#include "fast5_pool.hpp"
#include "chunk_pool.hpp"
#include "chunk.hpp"
#include "read_buffer.hpp"
#include "uncalled_opts.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(mapping, m) {
    m.doc() = "UNCALLED";

    py::class_<UncalledOpts>(m, "UncalledOpts")
        //Map constructor
        .def(py::init<const std::string &, //bwa_prefix,
                      const std::string &, //model_fname,
                      u32, // seed_len, 
                      u32, // min_aln_len,
                      u32, // min_rep_len, 
                      u32, // max_rep_copy, 
                      u32, // max_consec_stay,
                      u32, // max_paths, 
                      u32, // max_events_proc,
                      u32, // evt_winlen1,
                      u32, // evt_winlen2,
                      u16, // threads,
                      float, // evt_thresh1,
                      float, // evt_thresh2,
                      float, // evt_peak_height,
                      float, // evt_min_mean,
                      float, // evt_max_mean,
                      float, // max_stay_frac,
                      float, // min_seed_prob, 
                      float, // min_mean_conf,
                      float // min_top_conf);
                     >())


        //Simulate constructor
        .def(py::init<const std::string &, //bwa_prefix,
                      const std::string &, //model_fname,
                      u32, // seed_len, 
                      u32, // min_aln_len,
                      u32, // min_rep_len, 
                      u32, // max_rep_copy, 
                      u32, // max_consec_stay,
                      u32, // max_paths, 
                      u32, // max_events_proc,
                      u32, // max_chunks_proc,
                      u32, // evt_buffer_len,
                      u32, // evt_winlen1,
                      u32, // evt_winlen2,
                      u16, // threads,
                      u16, // num_channels,
                      u16, // chunk_len,
                      u16, // evt_batch_size,
                      float, // evt_timeout,
                      float, // evt_thresh1,
                      float, // evt_thresh2,
                      float, // evt_peak_height,
                      float, // evt_min_mean,
                      float, // evt_max_mean,
                      float, // max_stay_frac,
                      float, // min_seed_prob, 
                      float, // min_mean_conf,
                      float, // min_top_conf,
                      float // sim_speed);
                     >());

    //py::class_<Mapper>(m, "Mapper")
    //    .def(py::init<UncalledOpts &>())
    //    .def("map_fast5", &Mapper::map_fast5);

    py::class_<Paf> paf(m, "Paf");
    paf.def(py::init())
       .def("print", &Paf::print)
       .def("is_mapped", &Paf::is_mapped)
       .def("set_int", &Paf::set_int)
       .def("set_float", &Paf::set_float)
       .def("set_str", &Paf::set_str);

    py::enum_<Paf::Tag>(paf, "Tag")
        .value("MAP_TIME", Paf::Tag::MAP_TIME)
        .value("CHUNKS", Paf::Tag::CHUNKS)
        .value("UNBLOCK", Paf::Tag::UNBLOCK)
        .value("KEEP", Paf::Tag::KEEP)
        .export_values();

    //py::class_<Fast5Pool>(m, "Fast5Pool")
    //    .def(py::init<const UncalledOpts &>())
    //    .def("add_fast5s", &Fast5Pool::add_fast5s)
    //    .def("update", &Fast5Pool::update)
    //    .def("all_finished", &Fast5Pool::all_finished)
    //    .def("stop_all", &Fast5Pool::stop_all); 

    py::class_<ChunkPool>(m, "ChunkPool")
        .def(py::init<const UncalledOpts &>()) 
        .def("update", &ChunkPool::update)
        .def("all_finished", &ChunkPool::all_finished)
        .def("stop_all", &ChunkPool::stop_all)
        .def("add_chunk", &ChunkPool::add_chunk)
        .def("end_read", &ChunkPool::end_read);

    py::class_<Chunk>(m, "Chunk")
        .def(py::init<const std::string &, //id, 
                      u16, //channel
                      u32, //number, 
                      u64, //chunk_start, 
                      const std::string &, //dtype
                      const std::string & //raw_str
                     >())
        .def(py::init<const std::string&, //_id, 
                      u32, //number 
                      u16, //channel
                      u64, //chunk_start_sample, 
                      const std::vector<float> &, //raw_data, 
                      u32, //raw_st
                      u32  //raw_len
                     >())
        .def_static("set_calibration", &Chunk::set_calibration)
        .def("get_channel", &Chunk::get_channel)
        .def("get_number", &Chunk::get_number)
        .def("empty", &Chunk::empty)
        .def("print", &Chunk::print)
        .def("size", &Chunk::size);

    py::class_<Simulator>(m, "Simulator")
        .def(py::init<const UncalledOpts &>())
        .def("add_fast5s", &Simulator::add_fast5s)
        .def("get_read_chunks", &Simulator::get_read_chunks)
        .def("stop_receiving_read", &Simulator::stop_receiving_read)
        .def("unblock", &Simulator::unblock)
        .def("get_time", &Simulator::get_time)
        .def("start", &Simulator::start)
        .def("is_running", &Simulator::is_running);
    
    m.def("self_align", &self_align);
}

