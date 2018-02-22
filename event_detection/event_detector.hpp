#ifndef EVENT_DETECTOR_HPP
#    define EVENT_DETECTOR_HPP

extern "C" {
    #include "scrappie_structures.h"
}

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
} ;

class EventDetector {
    public:
    EventDetector(const detector_param &edparam);
    
    event_t add_signal(float s);

    private:
    size_t sig_center();
    float compute_tstat(size_t s, size_t w_length); 
    bool peak_detect(float current_value, Detector &detector);
    event_t create_event(); 

    const detector_param params;
    const size_t BUF_LEN;
    double *sum, *sumsq;

    size_t t, t_mod, evt_st;
    double evt_st_sum, evt_st_sumsq;

    Detector short_detector, long_detector;
};


//event_table detect_events(raw_table const rt, detector_param const edparam);

#endif                          /* EVENT_DETECTION_H */
