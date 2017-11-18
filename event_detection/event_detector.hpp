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

//class EventDetector {
//    EventDetector(const detector_param edparam);
//    
//    event_t add_signal(float s);
//};


event_table detect_events(raw_table const rt, detector_param const edparam);

#endif                          /* EVENT_DETECTION_H */
