#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "map_pool.hpp"
#include "self_align_ref.hpp"
#include "realtime_pool.hpp"
#include "client_sim.hpp"
#include "model_r94.inl"

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

    py::class_< PoreModel<KLEN> > model(m, "PoreModel");
    PoreModel<KLEN>::pybind_defs(model);
    m.attr("pmodel_r94_template") = py::cast(pmodel_r94_template);
    m.attr("pmodel_r94_complement") = py::cast(pmodel_r94_complement);

    m.def("self_align", &self_align);

    //BP operation functions
    m.def("kmer_count",    &kmer_count<KLEN>);
    m.def("str_to_kmer",   &str_to_kmer<KLEN>);
    m.def("kmer_comp",     &kmer_comp<KLEN>);
    m.def("kmer_revcomp",  &kmer_revcomp<KLEN>);
    m.def("kmer_head",     &kmer_head<KLEN>);
    m.def("kmer_base",     &kmer_base<KLEN>);
    m.def("kmer_to_str",   &kmer_to_str<KLEN>);
    m.def("seq_to_kmers",  &seq_to_kmers<KLEN>);
    m.def("kmer_neighbor", &kmer_neighbor<KLEN>);

}

