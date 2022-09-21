#ifndef _INCL_SIGNAL_PROCESSOR
#define _INCL_SIGNAL_PROCESSOR

#include <deque>
#include <valarray>
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
    i32 start, end;
    float scale, shift;
};

struct ProcessedRead {
    //std::vector<i32> event_starts, event_lengths;
    //std::vector<float> event_means, event_stdvs;
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


    std::pair<float,float> get_moments(size_t event_start, size_t event_end) {
        size_t n = event_end - event_start;
        std::valarray<float> event_means(n);

        size_t i = 0;
        for (size_t j = event_start; j < event_end; j++) {
            event_means[i++] = events[j].mean;
        }

        //auto mean = event_means.sum() / n;
        //auto deltas = event_means - mean;

        //const std::valarray<float>  filt_means = std::valarray<float>(event_means[std::abs(deltas) < 3.5*stdv]);
        //mean = filt_means.sum() / n;
        //auto deltas_f = filt_means - mean;
        //stdv = sqrt((deltas_f*deltas_f).sum() / n);

        std::sort(std::begin(event_means), std::end(event_means));
        auto median = event_means[event_means.size()/2];
        auto deltas = event_means - median;

        auto stdv = sqrt((deltas*deltas).sum() / n);

        //stdv = (event_means * event_means).sum() / event_means.size()

        //auto len = event_end - event_start;
        //float mean = 0, stdv = 0;
        //for (size_t i = event_start; i < event_end; i++) {
        //    mean += events[i].mean;
        //}
        //mean /= len;
        //
        ////for (auto &e : events) {
        //for (size_t i = event_start; i < event_end; i++) {
        //    float delta = events[i].mean - mean;
        //    stdv += delta*delta;
        //}
        //stdv = sqrt(stdv / len);

        //return {mean, stdv};
        return {median, stdv};
    }

    void normalize(NormVals prms) {
        assert(prms.start >= 0);
        assert(prms.end < events.size());
        for (size_t i = prms.start; i < prms.end; i++) {
            events[i].mean = events[i].mean * prms.scale + prms.shift;
        }
        norm.push_back(prms);
    }

    void normalize(float scale, float shift) {
        normalize({0, events.size(), scale, shift});
    }

    void normalize_mom(float tgt_mean, float tgt_stdv, size_t event_start, size_t event_end) {
        auto mom = get_moments(event_start, event_end);
        auto mean = mom.first, stdv = mom.second;

        NormVals norm;
        norm.start = event_start;
        norm.end = event_end;
        norm.scale = tgt_stdv / stdv;
        norm.shift = tgt_mean - norm.scale * mean;

        normalize(norm);
    }

    void normalize_mom(float tgt_mean, float tgt_stdv) {
        normalize_mom(tgt_mean, tgt_stdv, 0, events.size());
    }

    #ifdef PYBIND
    void set_events(py::array_t<Event> arr) {
        auto evts = PyArray<Event>(arr);
        events.clear();
        events.reserve(evts.size());
        for (size_t i = 0; i < evts.size(); i++) {
            events.push_back(evts[i]);
        }
    }
    #endif
};

template <typename ModelType>
class SignalProcessor {

    private:
    const ModelType &model_;
    EventDetector evdt_;
    Normalizer norm_;
    Normalizer::Params norm_prms_;

    public: 

    SignalProcessor(const ModelType &model, EventDetector::Params event_prms=EventDetector::PRMS_DEF, Normalizer::Params norm_prms=Normalizer::PRMS_DEF) : 
        model_(model),
        evdt_(event_prms),
        norm_prms_(norm_prms) {

        //set_norm_tgt(model.model_mean(), model.model_stdv());
    }

    ProcessedRead process(const ReadBuffer &read, bool normalize=true) {
        ProcessedRead ret = {};

        ret.events = evdt_.get_events(read.get_signal());
        
        if (normalize) {
            if (norm_prms_.mode == "model_mom") {
                ret.normalize_mom(norm_prms_.tgt_mean, norm_prms_.tgt_stdv);
            } else {
                throw std::runtime_error("Normalization mode not supported: " + norm_prms_.mode); 
            }
        }
        //for (auto &e : ret.events) {
        //    e.mean = e.mean * norm.scale + norm.shift;
        //}

        return ret;
    }

    void set_norm_tgt(float mean, float stdv) {
        norm_prms_.tgt_mean = mean;
        norm_prms_.tgt_stdv = stdv;
    }


    #ifdef PYBIND

    static void pybind(py::module_ &m, const std::string &suffix) {
        py::class_<SignalProcessor> s(m, ("SignalProcessor" + suffix).c_str());
        s.def(pybind11::init<const ModelType &, EventDetector::Params>());
        s.def("process", &SignalProcessor::process);
        s.def("set_norm_tgt", &SignalProcessor::set_norm_tgt);

    }
    #endif

};

#define PY_PROC_ARR(T, A, D) p.def_property_readonly(#A, \
    [](ProcessedRead &r) -> py::array_t<T> { \
        return py::array_t<T>(r.A.size(), r.A.data()); \
    }, D);


#ifdef PYBIND
void signal_processor_pybind(py::module_ &m) {

    PYBIND11_NUMPY_DTYPE(NormVals, start, end, scale, shift);

    //TODO move to ProcessedRead
    py::class_<ProcessedRead> p(m, "_ProcessedRead");
    p.def(pybind11::init());
    p.def(pybind11::init<const ProcessedRead &>());
    p.def("normalize", 
            static_cast<void (ProcessedRead::*)(float, float)> (&ProcessedRead::normalize),
            py::arg("scale"), py::arg("shift"));
    p.def("normalize_mom", 
            static_cast< void (ProcessedRead::*)(float, float)> (&ProcessedRead::normalize_mom),
            py::arg("tgt_mean"), py::arg("tgt_stdv"));
    p.def("normalize_mom", 
            static_cast< void (ProcessedRead::*)(float, float, size_t, size_t)> (&ProcessedRead::normalize_mom),
            py::arg("tgt_mean"), py::arg("tgt_stdv"), py::arg("start"), py::arg("end"));
    p.def_property_readonly("sample_start", &ProcessedRead::sample_start);
    p.def_property_readonly("sample_end", &ProcessedRead::sample_end);
    p.def("set_events", &ProcessedRead::set_events);
    PY_PROC_ARR(Event, events, "Un-normalized events");
    PY_PROC_ARR(NormVals, norm, "Normalizer values and read coordinates");
}
#endif

#endif
