/* 
 * This file was adptaed from https://github.com/nanoporetech/scrappie
 * Original resrion released under Mozzila Public Licence
 * Adapted for UNCALLED by Sam Kovaka <skovaka@gmail.com>
 */

#ifndef _INCL_EVENT_DETECTOR
#define _INCL_EVENT_DETECTOR

#include <vector>
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

    static const Params PRMS_DEF;

    Params PRMS;

    EventDetector(Params prms);
    EventDetector();

    ~EventDetector();
    
    void reset();
    bool add_sample(float s);
    Event get_event() const;
    std::vector<Event> get_events(const std::vector<float> &raw);

    float kmer_current() const;
    std::vector<float> get_means(const std::vector<float> &raw);

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
        PY_EVTD_METH(kmer_current, "");
        PY_EVTD_METH(get_means, "");
        PY_EVTD_METH(mean_event_len, "");

        evdt.def_readonly_static("PRMS_DEF", &EventDetector::PRMS_DEF, "");

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


#endif                          /* EVENT_DETECTION_H */
