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
    std::vector<Event> raw_events;
    std::vector<float> norm_events;
    std::vector<bool> mask;
    std::vector<NormVals> norm;
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
        ret.raw_events = evdt_.get_events(read.get_signal());
        ret.norm_events.reserve(ret.raw_events.size());
        for (auto &e : ret.raw_events) {
            ret.norm_events.push_back(e.mean);
        }
        norm_.set_signal(ret.norm_events);
        ret.norm.push_back({0, static_cast<i32>(ret.norm_events.size()), norm_.get_scale(), norm_.get_shift()});
        for (size_t i = 0; i < ret.norm_events.size(); i++) { 
            ret.norm_events[i] = norm_.pop();
        }

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
        PY_PROC_ARR(Event, raw_events, "Un-normalized events");
        PY_PROC_ARR(float, norm_events, "Normalized event means");
        PY_PROC_ARR(NormVals, norm, "Normalizer values and read coordinates");
    }
    #endif

};
#endif
