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

struct NormVals {
    u32 start, end;
    float scale, shift;
};

struct ProcessedRead {
    std::vector<Event> events;
    std::vector<NormVals> norm;
    //std::vector<bool> mask;

    u32 sample_start() const {
        return events[0].start;
    }

    u32 sample_end() const {
        auto e = events.back();
        return e.start + e.length;
    }

    void rescale(float scale, float shift) {
        normalize({0, sample_end(), scale, shift});
    }

    void normalize(NormVals prms) {
        //auto cmp = [](Event &e, u32 loc) -> bool {
        //    return e.start < loc;
        //};
        //auto st = std::lower_bound(events.begin(), events.end(), prms.start, cmp);
        //auto en = std::lower_bound(events.begin(), events.end(), prms.end, cmp);
        //for (auto e = st; e < en; e++) {
        norm = {prms};
        
        std::cout << "MILLIE QUIET " << prms.scale << " " << prms.shift << "\n";

        for (auto &e : events) {
            e.mean = e.mean * prms.scale + prms.shift;
        }
    }
};

template <typename ModelType>
class SignalProcessor {

    private:
    const ModelType &model_;
    EventDetector evdt_;
    Normalizer norm_;

    public: 

    SignalProcessor(const ModelType &model, EventDetector::Params event_prms=EventDetector::PRMS_DEF) : 
        model_(model),
        evdt_(event_prms) {
        norm_.set_target(model.model_mean(), model.model_stdv());
    }

    ProcessedRead process(const ReadBuffer &read) {
        ProcessedRead ret = {};

        ret.events = evdt_.get_events(read.get_signal());

        auto norm = norm_mom_params(ret.events);
        ret.normalize(norm);
        //for (auto &e : ret.events) {
        //    e.mean = e.mean * norm.scale + norm.shift;
        //}

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
        std::cout << "BLAH " 
                  << norm_.PRMS.tgt_mean << "\t"
                  << norm_.PRMS.tgt_stdv << "\t"
                  << mean << "\t"
                  << stdv << "\n";
        ret.scale = norm_.PRMS.tgt_stdv / stdv;
        ret.shift = norm_.PRMS.tgt_mean - ret.scale * mean;
        return ret;
    }

    #ifdef PYBIND

    #define PY_PROC_ARR(T, A, D) p.def_property_readonly(#A, \
        [](ProcessedRead &r) -> py::array_t<T> { \
            return py::array_t<T>(r.A.size(), r.A.data()); \
        }, D);

    static void pybind(py::module_ &m, const std::string &suffix) {
        py::class_<SignalProcessor> s(m, ("SignalProcessor" + suffix).c_str());
        s.def(pybind11::init<const ModelType &, EventDetector::Params>());
        s.def("process", &SignalProcessor::process);

    }
    #endif

};

#ifdef PYBIND
void signal_processor_pybind(py::module_ &m) {

    PYBIND11_NUMPY_DTYPE(NormVals, start, end, scale, shift);


    //TODO move to ProcessedRead
    py::class_<ProcessedRead> p(m, "_ProcessedRead");
    p.def(pybind11::init<const ProcessedRead &>());
    p.def("rescale", &ProcessedRead::rescale);
    p.def_property_readonly("sample_start", &ProcessedRead::sample_start);
    p.def_property_readonly("sample_end", &ProcessedRead::sample_end);
    PY_PROC_ARR(Event, events, "Un-normalized events");
    PY_PROC_ARR(NormVals, norm, "Normalizer values and read coordinates");

    SignalProcessor<PoreModel<5>>::pybind(m, "K5");
    SignalProcessor<PoreModel<10>>::pybind(m, "K10");
}
#endif

#endif
