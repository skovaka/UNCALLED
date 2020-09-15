#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "map_pool.hpp"
#include "self_align_ref.hpp"
#include "realtime_pool.hpp"
#include "client_sim.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_uncalled, m) {
    m.doc() = R"pbdoc(UNCALLED: a Utility for Nanopore Current ALignment to Large Expanses of DNA)pbdoc";

    py::class_<Conf> conf(m, "Conf");
    Conf::pybind_defs(conf);

    py::class_<MapPool> map_pool(m, "MapPool");
    MapPool::pybind_defs(map_pool);

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
        .value("DELAY", Paf::Tag::DELAY)
        .export_values();

    py::class_<BwaIndex<KLEN>> bwa_index(m, "BwaIndex");
    BwaIndex<KLEN>::pybind_defs(bwa_index);

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

    py::class_<Chunk> chunk(m, "Chunk");
    Chunk::pybind_defs(chunk);
    
    m.def("self_align", &self_align);

    py::class_<ReadBuffer>(m, "ReadBuffer")
        .def("empty",        &ReadBuffer::empty)
        .def("size",         &ReadBuffer::size)
        .def("get_id",       &ReadBuffer::get_id)
        .def("get_start",    &ReadBuffer::get_start)
        .def("get_end",      &ReadBuffer::get_end)
        .def("get_duration", &ReadBuffer::get_duration)
        .def("get_channel",  &ReadBuffer::get_channel)
        .def("get_raw",      &ReadBuffer::get_raw);

    py::class_<Fast5Reader>(m, "Fast5Reader")
        .def(py::init<u32, u32>())
        .def("add_fast5",       &Fast5Reader::add_fast5) 
        .def("load_fast5_list", &Fast5Reader::load_fast5_list)
        .def("add_read",        &Fast5Reader::add_read)
        .def("load_read_list",  &Fast5Reader::load_read_list)
        .def("pop_read",        &Fast5Reader::pop_read)
        .def("buffer_size",     &Fast5Reader::buffer_size)
        .def("fill_buffer",     &Fast5Reader::fill_buffer)
        .def("all_buffered",    &Fast5Reader::all_buffered)
        .def("empty",           &Fast5Reader::empty);

    py::class_<Event> event(m, "Event");
    py::class_<EventDetector> event_detector(m, "EventDetector");
    py::class_<EventDetector::Params> event_prms(event_detector, "Params");
    EventDetector::pybind_defs(event_detector, event_prms, event);

    py::class_<ClientSim> client_sim(m, "ClientSim");
    ClientSim::pybind_defs(client_sim);

    py::class_<Normalizer>(m, "Normalizer")
        .def(py::init<float, float>())
        .def("set_target", &Normalizer::set_target)
        .def("set_signal", &Normalizer::set_signal)
        .def("get_mean", &Normalizer::get_mean)
        .def("get_stdv", &Normalizer::get_stdv)
        .def("get_scale", &Normalizer::get_scale)
        .def("get_shift", &Normalizer::get_shift)
        .def("pop", &Normalizer::pop)
        .def("push", &Normalizer::push)
        .def("skip_unread", &Normalizer::skip_unread)
        .def("unread_size", &Normalizer::unread_size)
        .def("reset", &Normalizer::reset)
        .def("empty", &Normalizer::empty)
        .def("full", &Normalizer::full);
}

