#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ref_index.hpp"
#include "read_buffer_bc.hpp"
#include "signal_processor.hpp"
#include "dtw.hpp"
#include "intervals.hpp"
#include "eventalign.hpp"
#include "aln.hpp"

namespace py = pybind11;
using namespace pybind11::literals;


template<size_t K>
size_t pybind_kmer(py::module_ &m) {
    std::string suffix = "K"+std::to_string(K);
    using Model = PoreModel<K>;
    Model::pybind_defs(m, suffix);
    RefIndex<Model>::pybind_defs(m, suffix);//ref_index);
    BandedDTW<Model>::pybind_defs(m, suffix);
    //StaticBDTW<Model>::pybind_defs(m, suffix);
    GlobalDTW<Model>::pybind_defs(m, suffix);
    SignalProcessor<Model>::pybind(m, suffix);

    Sequence<Model>::pybind(m, suffix);
    Alignment<Model>::pybind(m, suffix);

    //m.def(("write_eventalign_"+suffix).c_str(), write_eventalign<PoreModel<K>>);
    auto fn = write_eventalign<Model>;
    m.def(("write_eventalign_"+suffix).c_str(), fn);
    return K;
}

template<size_t ...Ks>
std::vector<size_t> pybind_kmers(py::module_ &m) {
    return {(pybind_kmer<Ks>(m))...};
}

PYBIND11_MODULE(_uncalled, m) {
    m.doc() = R"pbdoc(UNCALLED: a Utility for Nanopore Current ALignment to Large Expanses of DNA)pbdoc";

    py::class_<Config> config(m, "_Conf");
    Config::pybind_defs(config);

    RefCoord::pybind_defs(m);

    py::class_<ReadBuffer> read_buffer(m, "ReadBuffer");
    ReadBuffer::pybind_defs(read_buffer);

    py::class_<ReadBufferBC, ReadBuffer> read_buffer_bc(m, "ReadBufferBC");
    ReadBufferBC::pybind_defs(read_buffer_bc);
    
    EventDetector::pybind_defs(m);

    py::class_<Normalizer> norm(m, "Normalizer");
    Normalizer::pybind_defs(norm);

    pybind_pore_model_params(m);

    py::bind_vector<std::vector<u8>>(m, "ArrayU8", py::buffer_protocol());
    py::bind_vector<std::vector<u16>>(m, "ArrayU16", py::buffer_protocol());
    py::bind_vector<std::vector<u32>>(m, "ArrayU32", py::buffer_protocol());
    pybind_kmers<5>(m);

    ProcessedRead::pybind(m);

    pybind_dtw(m);
    pybind_intervals(m);
    pybind_arrays(m);

    AlnDF::pybind(m);
    CmpDF::pybind(m);
}

