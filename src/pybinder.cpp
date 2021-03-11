#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "map_pool.hpp"
#include "self_align_ref.hpp"
#include "realtime_pool.hpp"
#include "map_pool_ord.hpp"
#include "client_sim.hpp"
#include "models.inl"
#include "dtw.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_uncalled, m) {
    m.doc() = R"pbdoc(UNCALLED: a Utility for Nanopore Current ALignment to Large Expanses of DNA)pbdoc";

    py::class_<Conf> conf(m, "Conf");
    Conf::pybind_defs(conf);

    py::class_<Mapper> mapper(m, "Mapper");
    Mapper::pybind_defs(mapper);

    py::class_<MapPool> map_pool(m, "MapPool");
    MapPool::pybind_defs(map_pool);

    py::class_<SeedTracker> seed_tracker(m, "SeedTracker");
    SeedTracker::pybind_defs(seed_tracker);

    py::class_<RealtimePool> realtime_pool(m, "RealtimePool");
    RealtimePool::pybind_defs(realtime_pool);

    py::class_<MapPoolOrd, RealtimePool> map_pool_ord(m, "MapPoolOrd");
    MapPoolOrd::pybind_defs(map_pool_ord);

    py::class_<ClientSim> client_sim(m, "ClientSim");
    ClientSim::pybind_defs(client_sim);

    py::class_<Paf> paf(m, "Paf");
    Paf::pybind_defs(paf);

    py::class_<BwaIndex<KLEN>> bwa_index(m, "BwaIndex");
    BwaIndex<KLEN>::pybind_defs(bwa_index);
    
    py::class_<ReadBuffer> read_buffer(m, "ReadBuffer");
    ReadBuffer::pybind_defs(read_buffer);
    
    py::class_<Fast5Read, ReadBuffer> fast5_read(m, "Fast5Read");
    Fast5Read::pybind_defs(fast5_read);

    py::class_<Fast5Reader> fast5_reader(m, "Fast5Reader");
    Fast5Reader::pybind_defs(fast5_reader);


    py::class_<Event> event(m, "Event");
    py::class_<EventDetector> event_detector(m, "EventDetector");
    EventDetector::pybind_defs(event_detector, event);

    py::class_<EventProfiler>  event_profiler(m, "EventProfiler");
    EventProfiler::pybind_defs(event_profiler);

    py::class_<Normalizer> norm(m, "Normalizer");
    Normalizer::pybind_defs(norm);

    py::class_< PoreModel<KLEN> > model(m, "PoreModel");
    PoreModel<KLEN>::pybind_defs(model);
    m.attr("PORE_MODELS") = py::cast(PORE_MODELS);
    //m.attr("pmodel_r94_dna_templ") = py::cast(pmodel_r94_dna_templ);
    //m.attr("pmodel_r94_dna_compl") = py::cast(pmodel_r94_dna_compl);
    //m.attr("pmodel_r94_rna_templ") = py::cast(pmodel_r94_rna_templ);
    //m.attr("pmodel_r94_rna_compl") = py::cast(pmodel_r94_rna_compl);

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

    //py::class_<DTW> dtw(m, "DTW");
    //DTW::pybind_defs(dtw);
    
    py::class_<DTWp> dtw_p(m, "DTWp");
    DTWp::pybind_defs(dtw_p);

    py::class_<DTWd> dtw_d(m, "DTWd");
    DTWd::pybind_defs(dtw_d);

    //py::class_<DTWr94p> dtw_r94p(m, "DTWr94p");
    //DTWr94p::pybind_defs(dtw_r94p);

    //py::class_<DTWr94d> dtw_r94d(m, "DTWr94d");
    //DTWr94d::pybind_defs(dtw_r94d);

    //py::class_<DTWr94pRNA> dtw_r94pRNA(m, "DTWr94pRNA");
    //DTWr94pRNA::pybind_defs(dtw_r94pRNA);

    //py::class_<DTWr94dRNA> dtw_r94dRNA(m, "DTWr94dRNA");
    //DTWr94dRNA::pybind_defs(dtw_r94dRNA);

    py::class_<DTWParams>dtwp(m, "DTWParams");
    dtwp.def_readwrite("dw", &DTWParams::dw);
    dtwp.def_readwrite("hw", &DTWParams::hw);
    dtwp.def_readwrite("vw", &DTWParams::vw);
    m.attr("DTW_EVENT_GLOB") = py::cast(DTW_EVENT_GLOB);
    m.attr("DTW_RAW_GLOB") = py::cast(DTW_RAW_GLOB);
    m.attr("DTW_GLOB") = py::cast(DTW_GLOB);
    m.attr("DTW_EVENT_QSUB") = py::cast(DTW_EVENT_QSUB);
    m.attr("DTW_EVENT_RSUB") = py::cast(DTW_EVENT_RSUB);
    m.attr("DTW_RAW_QSUB") = py::cast(DTW_EVENT_QSUB);
    m.attr("DTW_RAW_RSUB") = py::cast(DTW_EVENT_RSUB);
}

