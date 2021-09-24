#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "fast5_py_read.hpp"
#include "map_pool.hpp"
#include "self_align_ref.hpp"
#include "realtime_pool.hpp"
#include "map_pool_ord.hpp"
#include "simulator.hpp"
#include "dtw.hpp"
#include "dtw_banded.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

std::vector<bool> unpack_moves(u64 moves, u8 length) {
    std::vector<bool> ret(length);
    for (u32 i = 0; i < length; i++) {
        ret[i] = (moves >> i) & 1;
    }
    return ret;
}

PYBIND11_MODULE(_uncalled, m) {
    m.doc() = R"pbdoc(UNCALLED: a Utility for Nanopore Current ALignment to Large Expanses of DNA)pbdoc";

    py::class_<Config> config(m, "_Conf");
    Config::pybind_defs(config);

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

    py::class_<Simulator> simulator(m, "Simulator");
    Simulator::pybind_defs(simulator);

    py::class_<Paf> paf(m, "Paf");
    Paf::pybind_defs(paf);

    RefIndex<KLEN>::pybind_defs(m);//ref_index);
    RefCoord::pybind_defs(m);//ref_index);
    
    py::class_<ReadBuffer> read_buffer(m, "ReadBuffer");
    ReadBuffer::pybind_defs(read_buffer);

    py::class_<ReadBufferBC, ReadBuffer> read_buffer_bc(m, "ReadBufferBC");
    ReadBufferBC::pybind_defs(read_buffer_bc);
    
    py::class_<Fast5Read, ReadBufferBC> fast5_read(m, "Fast5Read");
    Fast5Read::pybind_defs(fast5_read);

    py::class_<Fast5PyRead, ReadBufferBC> fast5_py_read(m, "Fast5PyRead");
    Fast5PyRead::pybind_defs(fast5_py_read);

    py::class_<Fast5Reader> fast5_reader(m, "_Fast5Reader");
    Fast5Reader::pybind_defs(fast5_reader);

    py::class_<Fast5Iter, Fast5Reader> fast5_iter(m, "_Fast5Iter");
    Fast5Iter::pybind_defs(fast5_iter);

    py::class_<Fast5Dict, Fast5Reader> fast5_dict(m, "_Fast5Dict");
    Fast5Dict::pybind_defs(fast5_dict);

    py::class_<Event> event(m, "Event");
    py::class_<EventDetector> event_detector(m, "EventDetector");
    EventDetector::pybind_defs(event_detector, event);

    py::class_<EventProfiler>  event_profiler(m, "EventProfiler");
    EventProfiler::pybind_defs(event_profiler);

    py::class_<Normalizer> norm(m, "Normalizer");
    Normalizer::pybind_defs(norm);

    py::class_<Range> range(m, "Range");
    Range::pybind_defs(range);

    PoreModel<KLEN>::pybind_defs(m);

    m.def("self_align", &self_align);

    //BP operation functions
    auto nt_module = m.def_submodule("_nt", "Contains functions and constants for manipulating nucleotides and fixed length k-mers (5-mers) in binary format");
    nt_pybind_defs<KLEN>(nt_module);

    m.def("unpack_moves", &unpack_moves);

    //py::class_<DTW> dtw(m, "DTW");
    //DTW::pybind_defs(dtw);

    py::class_<BandedDTW> banded_dtw(m, "BandedDTW");
    BandedDTW::pybind_defs(banded_dtw);

    py::class_<StaticBDTW, BandedDTW> static_bdtw(m, "StaticBDTW");
    StaticBDTW::pybind_defs(static_bdtw);
    
    py::class_<DTWp> dtw_p(m, "DTWp");
    DTWp::pybind_defs(dtw_p);

    py::class_<DTWd, DTWp> dtw_d(m, "DTWd");
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

