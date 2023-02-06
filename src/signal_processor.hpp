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
        ProcessedRead ret = evdt_.process_read(read);//{norm_prms_};
		ret.norm_prms = norm_prms_;

        //ret.events = evdt_.get_events(read.get_signal());
		//ret.set_signal(read.get_signal());
        
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


#ifdef PYBIND
#endif

#endif
