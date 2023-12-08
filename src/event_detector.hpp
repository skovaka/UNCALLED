/* 
 * This file was adptaed from https://github.com/nanoporetech/scrappie
 * Original resrion released under Mozzila Public Licence
 * Adapted for UNCALLED by Sam Kovaka <skovaka@gmail.com>
 */

#ifndef _INCL_EVENT_DETECTOR
#define _INCL_EVENT_DETECTOR

#include <vector>
#include <unordered_map>
//#include "signal_processor.hpp"
#include "intervals.hpp"
#include "read_buffer.hpp"
#include "normalizer.hpp"
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

typedef struct {
    float mean; 
    float stdv; 
    u32 start;  
    u32 length; 
} Event;        

//
//#ifdef PYBIND
//PYBIND11_MAKE_OPAQUE(std::vector<Event>);
//#endif

struct NormVals {
    i32 start, end;
    float scale, shift;
};

struct ProcessedRead {
    //std::vector<i32> event_starts, event_lengths;
    //std::vector<float> event_means, event_stdvs;
    NormalizerParams norm_prms;
    std::vector<Event> events;
    std::vector<NormVals> norm;
    ValArray<float> signal;
    //std::vector<bool> mask;

    template <typename Container>
    void set_signal(Container signal_) {
        signal = ValArray<float>(&signal_[0], signal_.size());
    }

    size_t hard_mask(IntervalIndex<i32> sample_bounds) {
        std::vector<Event> masked;
        masked.reserve(events.size());
        size_t i = 0;
        for (size_t j = 0; j < sample_bounds.coords.size(); j++) {
        //for (auto &c : sample_bounds.coords) {
            auto &c = sample_bounds.coords[j];
            while (i < events.size() && (events[i].start + events[i].length) < c.start) i++;
            while (i < events.size() && events[i].start < c.end) {
                masked.push_back(events[i++]);
            }
        }
        events.swap(masked);
        return masked.size();
    }

    u32 sample_start() const {
        return events[0].start;
    }

    u32 sample_end() const {
        auto e = events.back();
        return e.start + e.length;
    }


    std::pair<float,float> get_moments(size_t event_start, size_t event_end) {
        size_t n = event_end - event_start;
        ValArray<float> event_means(n);

        size_t i = 0;
        for (size_t j = event_start; j < event_end; j++) {
            event_means[i++] = events[j].mean;
        }

        //const std::valarray<float>  filt_means = std::valarray<float>(event_means[std::abs(deltas) < 3.5*stdv]);

        float avg, dev;

        //if (norm_prms.median) {
        //    std::sort(std::begin(event_means), std::end(event_means));
        //    avg = event_means[n / 2];
        //    //mean = event_means.sum() / n;
        //} else {
            avg = event_means.sum() / n;
        //}

        auto deltas = event_means - avg;

        //if (norm_prms.median) {
        //    auto delta_abs = std::valarray<float>(std::abs(deltas));
        //    std::sort(std::begin(delta_abs), std::end(delta_abs));
        //    dev = delta_abs[n / 2];
        //} else {
            dev = sqrt((deltas*deltas).sum() / n);
        //}

        return {avg, dev};
    }

    void update_stdv(Event &e) const {
        auto slc = static_cast<std::valarray<float>>(signal[std::slice(e.start, e.length, 1)]);
        auto deltas = slc - e.mean;
        //e.mean = slc.sum() / e.length;
        e.stdv = sqrt((deltas*deltas).sum() / e.length);
        //e.mean = sqrt((deltas*deltas).sum() / e.length);
    }

    Event merge_events(size_t i, size_t j) const {
        if (i+1 >= j) {
            return events[i];
        }
        Event ret{0,0,events[i].start,0};
        for (auto k=i; k < j; k++) {
            ret.mean += events[k].mean*events[k].length;
            ret.length += events[k].length;
        }
        ret.mean /= ret.length;
        update_stdv(ret);
        return ret;
    }

    void normalize(NormVals prms) {
        assert(prms.start >= 0);
        assert(prms.end < events.size());

		signal = (signal * prms.scale) + prms.shift;

        //std::cout << "Normalizing " << prms.start << " " << prms.end << " " << prms.scale << " " << prms.shift << "\n";

        for (size_t i = prms.start; i < prms.end; i++) {
            auto &evt = events[i];
            evt.mean = evt.mean * prms.scale + prms.shift;
            update_stdv(evt);
        }
        //auto slc = std::slice(prms.start, prms.end-prms.start, 1);
		//signal[slc] = (static_cast<std::valarray<float>>(signal[slc]) * prms.scale) + prms.shift;
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

	#define PY_PROC_ARR(T, A, D) p.def_property_readonly(#A, \
    [](ProcessedRead &r) -> py::array_t<T> { \
        return py::array_t<T>(r.A.size(), r.A.data()); \
    }, D);

    #ifdef PYBIND
	static void pybind(py::module_ &m) {
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
		p.def("hard_mask", &ProcessedRead::hard_mask);
		p.def_readonly("signal", &ProcessedRead::signal);
		PY_PROC_ARR(Event, events, "Un-normalized events");        
		PY_PROC_ARR(NormVals, norm, "Normalizer values and read coordinates");    


	}

    void set_events(py::array_t<Event> arr) {
        auto evts = PyArray<Event>(arr);
        //events.clear();
        events.resize(evts.size());
        for (size_t i = 0; i < evts.size(); i++) {
            events[i] = evts[i];
        }
    }
    #endif
};

class EventDetector {

    public:
    typedef struct Params {
        u32 window_length1;
        u32 window_length2;
        float threshold1;
        float threshold2;
        float peak_height;
        float min_mean;
        float max_mean;
    } Params;

    static const Params PRMS_DEF, PRMS_450BPS, PRMS_70BPS;

    Params PRMS;

    EventDetector(Params prms);
    EventDetector();

    ~EventDetector();
    
    void reset();
    bool add_sample(float s);
    Event get_event() const;
    std::vector<Event> get_events(const ValArray<float> &raw);
    std::vector<Event> get_events2(const ValArray<float> &raw);

    ProcessedRead process_read(const ReadBuffer &read);
    ProcessedRead process_read_new(const ReadBuffer &read);

    float kmer_current() const;

    template <typename Container>
    std::vector<float> get_means(const Container &raw) {
        std::vector<float> events;
        events.reserve(raw.size() / PRMS.window_length2);
        reset();

        for (u32 i = 0; i < raw.size(); i++) {
            if (add_sample(raw[i])) {
                events.push_back(event_.mean);
            }
        }

        return events;
    }

    ValArray<float> get_means_py(py::array_t<float> raw);

    float mean_event_len() const;
    u32 event_to_bp(u32 evt_i, bool last=false) const;

    void set_calibration(float offset, float range, float digitisation);

    #ifdef DEBUG_EVDT
    struct Debug {
        std::vector<float> tstat1s, tstat2s;
        std::vector<u32> peak1_idxs, peak2_idxs;
        void reset() {
            tstat1s.clear();
            tstat2s.clear();
            peak1_idxs.clear();
            peak2_idxs.clear();
        }
    };
    Debug dbg_;

    Debug get_dbg() const {
        return dbg_;
    }
    #endif

    #ifdef PYBIND

    #ifdef PYDEBUG
    std::vector<Event> dbg_events_;
    #endif

    #define PY_EVTD_METH(P, D) evdt.def(#P, &EventDetector::P, D);
    #define PY_EVTD_PRM(P, D) prms.def_readwrite(#P, &EventDetector::Params::P, D);
    #define PY_EVT_VAL(P, D) evt.def_readwrite(#P, &Event::P, D);
    #define PY_DBG_VAL(P, D) dbg.def_readonly(#P, &Debug::P, D);

    static void pybind_defs(py::module_ &m) {
        py::class_<Event> evt(m, "Event");
        py::class_<EventDetector> evdt(m, "EventDetector");

        PYBIND11_NUMPY_DTYPE(Event, start, length, mean, stdv);

        evdt.def(pybind11::init<Params>());
        evdt.def(pybind11::init());
        PY_EVTD_METH(reset, "");
        PY_EVTD_METH(add_sample, "");
        //:PY_EVTD_METH(get_event, "");
        PY_EVTD_METH(get_events, "");
        PY_EVTD_METH(process_read, "");
        PY_EVTD_METH(process_read_new, "");
        PY_EVTD_METH(kmer_current, "");
        
        evdt.def("get_means", &EventDetector::get_means_py);
        //PY_EVTD_METH(get_means, "");
        PY_EVTD_METH(mean_event_len, "");

        evdt.def_readonly_static("PRMS_DEF", &EventDetector::PRMS_DEF, "");
        evdt.def_readonly_static("PRMS_450BPS", &EventDetector::PRMS_450BPS, "");
        evdt.def_readonly_static("PRMS_70BPS", &EventDetector::PRMS_70BPS, "");

        pybind11::class_<Params> prms(evdt, "Params", "");
        prms.def(pybind11::init(), "");
        prms.def(pybind11::init<Params>(), "");
        PY_EVTD_PRM(window_length1, "EventDetector t-stat window length 1");
        PY_EVTD_PRM(window_length2, "EventDetector t-stat window length 2");
        PY_EVTD_PRM(threshold1, "EventDetector t-stat threshold 1");
        PY_EVTD_PRM(threshold2, "EventDetector t-stat threshold 2");
        PY_EVTD_PRM(peak_height, "EventDetector peak_height");
        PY_EVTD_PRM(min_mean, "Minimum mean event pA");
        PY_EVTD_PRM(max_mean, "Maximum mean event pA");

        //TODO define as numpy records, pass vectors of events
        PY_EVT_VAL(mean, "");
        PY_EVT_VAL(stdv, "");
        PY_EVT_VAL(start, "");
        PY_EVT_VAL(length, "");

        #ifdef DEBUG_EVDT
        evdt.def_readonly("dbg", &EventDetector::dbg_);
        pybind11::class_<Debug> dbg(evdt, "Debug");
        PY_DBG_VAL(tstat1s)
        PY_DBG_VAL(tstat2s)
        PY_DBG_VAL(peak1_idxs)
        PY_DBG_VAL(peak2_idxs)
        #endif

    }

    #endif

    private:

    struct Detector {
        i32 DEF_PEAK_POS;
        float DEF_PEAK_VAL;
        float threshold;
        u32 window_length;
        u32 masked_to;
        i32 peak_pos;
        float peak_value;
        bool valid_peak;
    };

    u32 get_buf_mid();
    float compute_tstat(u32 w_length); 
    bool peak_detect(float current_value, Detector &detector);
    Event create_event(u32 evt_en); 
    float calibrate(float v);

    const u32 BUF_LEN;
    double *sum, *sumsq;

    u32 t, buf_mid, evt_st;
    double evt_st_sum, evt_st_sumsq;

    float cal_offset_, cal_coef_;

    Event event_;
    float len_sum_;
    u32 total_events_;

    Detector short_detector, long_detector;
};


//class EventDetectorMeta : public EventDetector {
//    struct Dbg {
//        std::unordered_map<int, std::vector<float>> tstats;
//        std::unordered_map<int, std::vector<u32>> peaks;
//
//        void init_window(u32 window_length) {
//            tstats[window_length] = {};
//            peaks[window_length] = {};
//        }
//
//        void reset() {
//            for (auto &i : tstats) i.second.clear();
//            for (auto &i : peaks) i.second.clear();
//        }
//    };
//
//    EventDetectorMeta();
//    void reset();
//    float compute_tstat(u32 w_length);
//    peak_detect(float current_value, EventDetector::Detector &detector);
//};


#endif                          /* EVENT_DETECTION_H */
