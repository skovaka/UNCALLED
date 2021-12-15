#ifndef _INCL_SIGNAL_PROCESSOR
#define _INCL_SIGNAL_PROCESSOR

#include <deque>
#include "read_buffer.hpp"
#include "event_detector.hpp"
#include "normalizer.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

const KmerLen SIG_KLEN = KmerLen::k5;

struct NormVals {
    i32 start, end;
    float scale, shift;
};

struct ProcessedRead {
    std::vector<Event> events;
    std::vector<NormVals> norm;
    //std::vector<bool> mask;
};

class SignalProcessor {

    private:
    const PoreModel<SIG_KLEN> &model_;
    EventDetector evdt_;
    Normalizer norm_;

    public: 

    SignalProcessor(const PoreModel<SIG_KLEN> &model, EventDetector::Params event_prms=EventDetector::PRMS_DEF) : 
        model_(model),
        evdt_(event_prms) {
        norm_.set_target(model.model_mean(), model.model_stdv());
    }

    ProcessedRead process(const ReadBuffer &read) {
        ProcessedRead ret = {};

        ret.events = evdt_.get_events(read.get_signal());

        auto norm = norm_mom_params(ret.events);
        for (auto &e : ret.events) {
            e.mean = e.mean * norm.scale + norm.shift;
        }

        return ret;
    }

    NormVals norm_mom_params(const std::vector<Event> &events) {
        float mean = 0, stdv = 0;
        for (auto &e : events) {
            mean += e.mean;
        }
        mean /= events.size();
        
        for (auto &e : events) {
            float delta = e.mean - mean;
            stdv += delta*delta;
        }
        stdv = sqrt(stdv / events.size());

        NormVals ret;
        ret.start = 0;
        ret.end = events.size();
        ret.scale = norm_.PRMS.tgt_stdv / stdv;
        ret.shift = norm_.PRMS.tgt_mean - ret.scale * mean;
        return ret;
    }

    #ifdef PYBIND

    #define PY_PROC_ARR(T, A, D) p.def_property_readonly(#A, \
        [](ProcessedRead &r) -> py::array_t<T> { \
            return py::array_t<T>(r.A.size(), r.A.data()); \
        }, D);

    static void pybind_defs(py::module_ &m) {
        py::class_<SignalProcessor> s(m, "_SignalProcessor");
        s.def(pybind11::init<const PoreModel<SIG_KLEN> &, EventDetector::Params>());
        s.def("process", &SignalProcessor::process);

        PYBIND11_NUMPY_DTYPE(NormVals, start, end, scale, shift);

        py::class_<ProcessedRead> p(m, "_ProcessedRead");
        p.def(pybind11::init<const ProcessedRead &>());
        PY_PROC_ARR(Event, events, "Un-normalized events");
        PY_PROC_ARR(NormVals, norm, "Normalizer values and read coordinates");
    }
    #endif

};
#endif
