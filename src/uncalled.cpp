#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "self_align_ref.hpp"
#include "mapper.hpp"
//#include "simulator.hpp"
#include "map_pool.hpp"
#include "realtime_pool.hpp"
#include "chunk.hpp"
#include "read_buffer.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(mapping, m) {
    m.doc() = "UNCALLED";

    py::class_<Conf> conf(m, "Conf");

    conf.def(py::init<const std::string &>())
        .def("load_conf", &Conf::load_conf);
    Conf::add_pybind_vars(conf);

    py::class_<Paf> paf(m, "Paf");
    paf.def(py::init())
       .def("print_paf", &Paf::print_paf)
       .def("is_mapped", &Paf::is_mapped)
       .def("is_ended", &Paf::is_ended)
       .def("set_int", &Paf::set_int)
       .def("set_float", &Paf::set_float)
       .def("set_str", &Paf::set_str);

    py::enum_<Paf::Tag>(paf, "Tag")
        .value("MAP_TIME", Paf::Tag::MAP_TIME)
        .value("EJECT", Paf::Tag::EJECT)
        .value("IN_SCAN", Paf::Tag::IN_SCAN)
        .value("ENDED", Paf::Tag::ENDED)
        .value("KEEP", Paf::Tag::KEEP)
        .export_values();

    //.def(py::init<const std::string &,const std::string &,u32>())

    py::class_<MapPool>(m, "MapPool")
        .def(py::init<Conf &>())
        .def("update", &MapPool::update)
        .def("add_fast5", &MapPool::add_fast5)
        .def("running", &MapPool::running)
        .def("stop", &MapPool::stop); 

    py::class_<RealtimePool>(m, "RealtimePool")
        .def(py::init<Conf &>()) 
        .def("update", &RealtimePool::update)
        .def("all_finished", &RealtimePool::all_finished)
        .def("stop_all", &RealtimePool::stop_all)
        .def("add_chunk", &RealtimePool::add_chunk);
    
    py::enum_<RealtimeParams::Mode>(m, "RealtimeMode")
        .value("DEPLETE", RealtimeParams::Mode::DEPLETE)
        .value("ENRICH", RealtimeParams::Mode::ENRICH)
        .export_values();

    py::enum_<RealtimeParams::ActiveChs>(m, "ActiveChs")
        .value("FULL", RealtimeParams::ActiveChs::FULL)
        .value("EVEN", RealtimeParams::ActiveChs::EVEN)
        .value("ODD", RealtimeParams::ActiveChs::ODD)
        .export_values();

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
        .def("get_channel", &Chunk::get_channel)
        .def("get_number", &Chunk::get_number)
        .def("empty", &Chunk::empty)
        .def("print", &Chunk::print)
        .def("size", &Chunk::size);

    //py::class_<Simulator>(m, "Simulator")
    //    .def(py::init())
    //    .def("add_fast5s", &Simulator::add_fast5s)
    //    .def("get_read_chunks", &Simulator::get_read_chunks)
    //    .def("stop_receiving_read", &Simulator::stop_receiving_read)
    //    .def("unblock", &Simulator::unblock)
    //    .def("get_time", &Simulator::get_time)
    //    .def("start", &Simulator::start)
    //    .def("is_running", &Simulator::is_running);
    
    m.def("self_align", &self_align);
}

