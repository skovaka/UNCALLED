#ifdef PYBIND

#ifndef ALN_HPP
#define ALN_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "pore_model.hpp"
namespace py = pybind11;

struct AlnDF {
    IntervalIndex<i64> index;
    IntervalIndex<i32> samples;
    ValArray<float> current, current_sd; 

    AlnDF() {}

    AlnDF(IntervalIndex<i64> index_) : index(index_) {
        //current(index.length),
        //current_sd(index.length)
        current = ValArray<float>(index_.length);
        current_sd = ValArray<float>(index_.length);
    }

    AlnDF(IntervalIndex<i64> index_, py::array_t<i32> start_, py::array_t<i32> length_, py::array_t<float> current_, py::array_t<float> current_sd_) : 
        index(index_),
        samples(start_, length_),
        current(init_arr(current_)),
        current_sd(init_arr(current_sd_)) {}

    AlnDF(IntervalIndex<i64> index_, IntervalIndex<i32> &samples_, py::array_t<float> current_, py::array_t<float> current_sd_) : 
        index(index_),
        samples(samples_),
        current(init_arr(current_)),
        current_sd(init_arr(current_sd_)) {}

    AlnDF slice(size_t i, size_t j) {
        AlnDF ret(index.islice(i, j));
        for (size_t k = i; k < j; k++) {
            ret.samples.append(samples.coords[k]);
            ret.current[k-i] = current[k];
            ret.current_sd[k-i] = current_sd[k];
        }
        return ret;
    }

    bool empty() const {
        return size() == 0;
    }

    size_t size() const {
        return index.length;
    }

    template <typename T>
    static ValArray<T> init_arr(py::array_t<T> &a) {
        auto info = a.request();
        return ValArray<T>(static_cast<T*>(info.ptr), static_cast<size_t>(info.shape[0]));
    }

    static py::class_<AlnDF> pybind(py::module_ &m) {
        py::class_<AlnDF> c(m, "_AlnDF");
        c.def(py::init<IntervalIndex<i64> &>());
        c.def(py::init<IntervalIndex<i64>&, IntervalIndex<i32>&, py::array_t<float>, py::array_t<float>>());
        c.def(py::init<IntervalIndex<i64>&, py::array_t<i32>, py::array_t<i32>, py::array_t<float>, py::array_t<float>>());
        c.def("slice", &AlnDF::slice);
        c.def("__len__", &AlnDF::size);
        c.def_readwrite("index", &AlnDF::index);
        c.def_readwrite("samples", &AlnDF::samples);
        c.def_readwrite("current", &AlnDF::current);
        c.def_readwrite("current_sd", &AlnDF::current_sd);
        return c;
    }
};

template <typename ModelType>
struct Alignment {
    std::string read_id;
    Sequence<ModelType> seq;
    AlnDF dtw, moves;

    Alignment(const std::string &read_id_, Sequence<ModelType> seq_) :
        read_id(read_id_), seq(seq_) {
    }

    void set_dtw(AlnDF df) {
        dtw = df;
    }

    void set_moves(AlnDF df) {
        moves = df;
    }

    static void pybind(py::module_ &m, std::string suffix) {
        py::class_<Alignment> c(m, ("_Alignment"+suffix).c_str());
        c.def(py::init<const std::string &, Sequence<ModelType>>());
        c.def("set_dtw", &Alignment::set_dtw);
        c.def("set_moves", &Alignment::set_moves);
        c.def_readonly("read_id", &Alignment::read_id);
        c.def_readonly("seq", &Alignment::seq);
        c.def_readonly("_dtw", &Alignment::dtw);
        c.def_readonly("_moves", &Alignment::moves);
    }
};


#endif
#endif
