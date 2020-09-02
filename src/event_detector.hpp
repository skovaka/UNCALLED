/* 
 * This file was adptaed from https://github.com/nanoporetech/scrappie
 * Original resrion released under Mozzila Public Licence
 * Adapted for UNCALLED by Sam Kovaka <skovaka@gmail.com>
 */

#ifndef _INCL_EVENT_DETECTOR
#define _INCL_EVENT_DETECTOR

#include <fast5.hpp>
#include "util.hpp"

//typedef fast5::EventDetection_Event Event;
typedef fast5::Basecall_Event BCEvent;
typedef fast5::Raw_Sample RawSample;

typedef struct {
    float mean;
    float stdv;
    u32 start;
    u32 length;
} Event;

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
    bool add_sample(RawSample s);
    Event get_event() const;
    std::vector<Event> get_events(const std::vector<RawSample> &raw);

    float get_mean() const;
    std::vector<float> get_means(const std::vector<RawSample> &raw);

    float mean_event_len() const;
    u32 event_to_bp(u32 evt_i, bool last=false) const;

    private:
    u32 get_buf_mid();
    float compute_tstat(u32 w_length); 
    bool peak_detect(float current_value, Detector &detector);
    Event create_event(u32 evt_en); 

    const u32 BUF_LEN;
    double *sum, *sumsq;

    u32 t, buf_mid, evt_st;
    double evt_st_sum, evt_st_sumsq;

    Event event_;
    float len_sum_;
    u32 total_events_;

    Detector short_detector, long_detector;
};


#endif                          /* EVENT_DETECTION_H */
