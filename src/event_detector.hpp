/* 
 * This file was adptaed from https://github.com/nanoporetech/scrappie
 * Original resrion released under Mozzila Public Licence
 * Adapted for UNCALLED by Sam Kovaka <skovaka@gmail.com>
 */

#ifndef _INCL_EVENT_DETECTOR
#define _INCL_EVENT_DETECTOR

#include <fast5.hpp>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#endif


typedef struct {
    float mean;
    float stdv;
    u32 start;
    u32 length;
} Event;

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

    static Params const PRMS_DEF;

    Params PRMS;

    EventDetector(Params prms);
    EventDetector();

    ~EventDetector();
    
    void reset();
    bool add_sample(float s);
    Event get_event() const;
    std::vector<Event> get_events(const std::vector<float> &raw);

    float get_mean() const;
    std::vector<float> get_means(const std::vector<float> &raw);

    float mean_event_len() const;
    u32 event_to_bp(u32 evt_i, bool last=false) const;

    void set_calibration(float offset, float range, float digitisation);

    #ifdef PYBIND

    #define PY_EVTD_METH(P) d.def(#P, &EventDetector::P);
    #define PY_EVTD_PRM(P) p.def_readwrite(#P, &EventDetector::Params::P);
    #define PY_EVT_VAL(P) e.def_readwrite(#P, &Event::P);

    static void pybind_defs(
            pybind11::class_<EventDetector> &d,
            pybind11::class_<Event> &e) {

        d.def(pybind11::init<Params>());
        d.def(pybind11::init());
        PY_EVTD_METH(reset);
        PY_EVTD_METH(add_sample);
        PY_EVTD_METH(get_event);
        PY_EVTD_METH(get_events);
        PY_EVTD_METH(get_mean);
        PY_EVTD_METH(get_means);
        PY_EVTD_METH(mean_event_len);

        d.def_readonly_static("PRMS_DEF", &EventDetector::PRMS_DEF);

        pybind11::class_<Params> p(d, "Params");
        p.def(pybind11::init());
        p.def(pybind11::init<Params>());
        PY_EVTD_PRM(window_length1);
        PY_EVTD_PRM(window_length2);
        PY_EVTD_PRM(threshold1);
        PY_EVTD_PRM(threshold2);
        PY_EVTD_PRM(peak_height);
        PY_EVTD_PRM(min_mean);
        PY_EVTD_PRM(max_mean);

        PY_EVT_VAL(mean);
        PY_EVT_VAL(stdv);
        PY_EVT_VAL(start);
        PY_EVT_VAL(length);
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
