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

    py::class_<RealtimePool> realtime_pool(m, "RealtimePool");
    RealtimePool::pybind_defs(realtime_pool);

    py::class_<ClientSim> client_sim(m, "ClientSim");
    ClientSim::pybind_defs(client_sim);

    py::class_<Paf> paf(m, "Paf");
    Paf::pybind_defs(paf);

    py::class_<BwaIndex<KLEN>> bwa_index(m, "BwaIndex");
    BwaIndex<KLEN>::pybind_defs(bwa_index);

    py::class_<Chunk> chunk(m, "Chunk");
    Chunk::pybind_defs(chunk);
    
    py::class_<ReadBuffer> read_buffer(m, "ReadBuffer");
    ReadBuffer::pybind_defs(read_buffer);

    py::class_<Fast5Reader> fast5_reader(m, "Fast5Reader");
    Fast5Reader::pybind_defs(fast5_reader);

    py::class_<Event> event(m, "Event");
    py::class_<EventDetector> event_detector(m, "EventDetector");
    EventDetector::pybind_defs(event_detector, event);

    py::class_<Normalizer> norm(m, "Normalizer");
    Normalizer::pybind_defs(norm);

    m.def("self_align", &self_align);
}

