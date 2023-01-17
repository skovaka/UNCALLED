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



template <typename ModelType>
class SignalProcessor {

    private:
    const ModelType &model_;
    EventDetector evdt_;
    Normalizer norm_;
    NormalizerParams norm_prms_;

    public: 

    SignalProcessor(const ModelType &model, EventDetector::Params event_prms=EventDetector::PRMS_DEF, NormalizerParams norm_prms=NORMALIZER_PRMS_DEF) : 
        model_(model),
        evdt_(event_prms),
        norm_prms_(norm_prms) {

        //set_norm_tgt(model.model_mean(), model.model_stdv());
    }

    ProcessedRead process(const ReadBuffer &read, bool normalize=true) {
        ProcessedRead ret = {norm_prms_};

        ret.events = evdt_.get_events(read.get_signal());
        
        if (normalize) {
            if (norm_prms_.mode == "model_mom" || norm_prms_.mode == "ref_mom") {
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
        s.def(pybind11::init<const ModelType &, EventDetector::Params, NormalizerParams>());
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
