#ifndef EVENT_DETECTOR_HPP
#define EVENT_DETECTOR_HPP

#include "fast5.hpp"

typedef fast5::EventDetection_Event Event;
typedef fast5::Basecall_Event BCEvent;
typedef fast5::Raw_Sample RawSample;

//extern "C" {
//    #include "scrappie_structures.h"
//}
//

typedef struct {
    size_t window_length1;
    size_t window_length2;
    float  threshold1;
    float  threshold2;
    float  peak_height;
} detector_param;

static detector_param const event_detection_defaults = {
    .window_length1 = 3,
    .window_length2 = 6,
    .threshold1     = 1.4f,
    .threshold2     = 9.0f,
    .peak_height    = 0.2f
};

struct Detector {
    int DEF_PEAK_POS;
    float DEF_PEAK_VAL;
    float threshold;
    size_t window_length;
    size_t masked_to;
    int peak_pos;
    float peak_value;
    bool valid_peak;
};

class EventDetector {
    public:
    EventDetector(const detector_param &edparam, double min_mean, double max_mean);
    ~EventDetector();
    
    void reset();
    Event add_sample(RawSample s);
    std::vector<Event> get_all_events(std::vector<RawSample> raw);

    private:
    size_t get_buf_mid();
    float compute_tstat(size_t w_length); 
    bool peak_detect(float current_value, Detector &detector);
    Event create_event(int evt_en); 

    const detector_param params;
    const size_t BUF_LEN;
    const double MIN_MEAN, MAX_MEAN;
    double *sum, *sumsq;

    size_t t, buf_mid, evt_st;
    double evt_st_sum, evt_st_sumsq;

    Detector short_detector, long_detector;
};


#endif                          /* EVENT_DETECTION_H */
