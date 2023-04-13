/* 
 * This file was adptaed from https://github.com/nanoporetech/scrappie
 * Original vesrion released under Mozzila Public Licence
 * Adapted for UNCALLED by Sam Kovaka <skovaka@gmail.com>
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "event_detector.hpp"

const EventDetector::Params 
    EventDetector::PRMS_DEF = {
        window_length1 : 0,
        window_length2 : 0,
        threshold1     : -1,
        threshold2     : -1,
        peak_height    : -1,
        min_mean       : -1,
        max_mean       : -1
    }, EventDetector::PRMS_450BPS = {
        window_length1 : 3,
        window_length2 : 6,
        threshold1     : 1.4,
        threshold2     : 9.0,
        peak_height    : 0.2,
        min_mean       : -200,
        max_mean       : 200
    }, EventDetector::PRMS_70BPS = {
        window_length1 : 7,
        window_length2 : 12,
        threshold1     : 2.8,
        threshold2     : 18.0,
        peak_height    : 0.2,
        min_mean       : -200,
        max_mean       : 200
    };

EventDetector::Params init_params(EventDetector::Params prms, const EventDetector::Params &defs=EventDetector::PRMS_450BPS) {
    #define SET_DEFAULT(N) if (prms.N == EventDetector::PRMS_DEF.N) prms.N = defs.N;
    SET_DEFAULT(window_length1)
    SET_DEFAULT(window_length2)
    SET_DEFAULT(threshold1)
    SET_DEFAULT(threshold2)
    SET_DEFAULT(peak_height)
    SET_DEFAULT(min_mean)
    SET_DEFAULT(max_mean)
    return prms;
}

EventDetector::EventDetector(Params prms) :
    PRMS(init_params(prms)),
    BUF_LEN (1 + PRMS.window_length2 * 2),
    cal_offset_(0),
    cal_coef_(1) {

    sum = new double[BUF_LEN];
    sumsq = new double[BUF_LEN];

    reset();
}

EventDetector::EventDetector() : EventDetector(PRMS_DEF) {}

EventDetector::~EventDetector() {
    delete[] sum;
    delete[] sumsq;
}

void EventDetector::reset() {
    //for (size_t i = 0; i < BUF_LEN; i++) {
    //    sum[i] = sumsq[i] = 0;
    //}
    sum[0] = sumsq[0] = 0.0;
    t = 1;
    evt_st = 0;
    evt_st_sum = evt_st_sumsq = 0.0;

    len_sum_ = 0;
    total_events_ = 0;

    //TODO store defaults elsewhere
    short_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .threshold = PRMS.threshold1,
        .window_length = PRMS.window_length1,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false
    };

    long_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .threshold = PRMS.threshold2,
        .window_length = PRMS.window_length2,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false
    };

    #ifdef DEBUG_EVDT
    dbg_.reset();
    #endif
}

u32 EventDetector::get_buf_mid() {
    return t - (BUF_LEN / 2) - 1;
}

bool EventDetector::add_sample(float s) {

    u32 t_mod = t % BUF_LEN;
    
    if (t_mod > 0) {
        sum[t_mod] = sum[t_mod-1] + s;
        sumsq[t_mod] = sumsq[t_mod-1] + s*s;
    } else {
        sum[t_mod] = sum[BUF_LEN-1] + s;
        sumsq[t_mod] = sumsq[BUF_LEN-1] + s*s;
    }

    t++;
    buf_mid = get_buf_mid();

    double tstat1 = compute_tstat(PRMS.window_length1),
           tstat2 = compute_tstat(PRMS.window_length2);

    //std::cout << tstat1 << "\t" << tstat2 << "\n";

    bool p1 = peak_detect(tstat1, short_detector),
         p2 = peak_detect(tstat2, long_detector);

    #ifdef DEBUG_EVDT
    if (buf_mid < t) {
        if (p1) dbg_.peak1_idxs.push_back(buf_mid);
        if (p2) dbg_.peak2_idxs.push_back(buf_mid);
        dbg_.tstat1s.push_back(tstat1);
        dbg_.tstat2s.push_back(tstat2);
    }
    #endif

    if (p1 || p2) {
        create_event(buf_mid-(PRMS.window_length1/2+1)); //+3?

        if (event_.mean >= PRMS.min_mean && event_.mean <= PRMS.max_mean) {
            return true;
        }
    }
    
    return false;
}

std::vector<Event> EventDetector::get_events(const ValArray<float> &raw) {
    std::vector<Event> events;
    events.reserve(raw.size() / PRMS.window_length2);
    reset();

    for (u32 i = 0; i < raw.size(); i++) {
        if (add_sample(raw[i])) {
            events.push_back(event_);
        }
    }

    return events;
}

std::vector<Event> EventDetector::get_events2(const ValArray<float> &raw) {
    std::vector<Event> events;
    events.reserve(raw.size() / PRMS.window_length2);
    reset();

    float mean = 0, stdv = 0;

    for (u32 i = 0; i < raw.size(); i++) {
        if (add_sample(raw[i])) {
            events.push_back(event_);
            mean += event_.mean;
        }
    }
    
    mean /= events.size();

    for (auto &e : events) {
        auto delta = e.mean - mean;
        stdv += delta*delta;
    }

    stdv = sqrt(stdv / events.size());

    auto win = stdv * 2,
         min_mean = mean - win, 
         max_mean = mean + win;

    //std::cout << min_mean << " " << max_mean << " BLAH\n";
    size_t i = 0, j = 0;
    for (; i < events.size(); i++) {
        if (events[i].mean >= min_mean && events[i].mean <= max_mean) {
            events[j++] = events[i];
        }
    }
    events.resize(j);

    //IntervalArray<i32> coords;
    //coords.reserve(events.size();
    //for (auto &e : events) {
    //    if (e.mean >= min_mean && e.mean <= max_mean) {
    //        coords.append(e.start, e.start+e.length);
    //    }
    //}

    return events;
}

ProcessedRead EventDetector::process_read(const ReadBuffer &read) {
    ProcessedRead ret;
    ret.events = get_events(read.signal);
    ret.set_signal(read.signal);
    return ret;
}

ProcessedRead EventDetector::process_read_new(const ReadBuffer &read) {
    ProcessedRead ret;
    ret.events = get_events2(read.signal);
    ret.set_signal(read.signal);
    return ret;
}

Event EventDetector::get_event() const {
    return event_;
}

ValArray<float> EventDetector::get_means_py(py::array_t<float> raw_py) {
    PyArray<float> raw(raw_py);
    auto events = get_means(raw);
    ValArray<float> ret(events.data(), events.size());
    return ret;
}

//TODO: template with float, double, Event?
//template <typename Container>
//std::vector<float> EventDetector::get_means(const Container &raw) {
//}

float EventDetector::kmer_current() const {
    return event_.mean;
}

float EventDetector::mean_event_len() const {
    return len_sum_ / total_events_;
}

void EventDetector::set_calibration(float offset, float range, float digitisation) {
    cal_offset_ = offset;
    cal_coef_ = range / digitisation;
}

float EventDetector::calibrate(float v) {
    return (v + cal_offset_) * cal_coef_;
}

/**
 *   Compute windowed t-statistic from summary information
 *
 *   @param sum       double[d_length]  Cumulative sums of data (in)
 *   @param sumsq     double[d_length]  Cumulative sum of squares of data (in)
 *   @param d_length                    Length of data vector
 *   @param w_length                    Window length to calculate t-statistic over
 *
 *   @returns float array containing tstats.  Returns NULL on error
 **/
float EventDetector::compute_tstat(u32 w_length) {
    assert(w_length > 0);

    //float *tstat = (float *) calloc(d_length, sizeof(float));

    const float eta = FLT_MIN;
    const float w_lengthf = (float) w_length;

    // Quick return:
    //   t-test not defined for number of points less than 2
    //   need at least as many points as twice the window length
    if (t < BUF_LEN || w_length < 2) {
    //if (t < w_length*2 || w_length < 2) {
        return 0;
    }

    // fudge boundaries
    //for (u32 i = 0; i < w_length; ++i) {
    //    tstat[i] = 0;
    //    tstat[d_length - i - 1] = 0;
    //}

    u32 i = buf_mid % BUF_LEN,
           st = (buf_mid - w_length) % BUF_LEN,
           en = (buf_mid + w_length) % BUF_LEN;

    //std::cout << i << " " << st << " " << en << "\n";

    //double sum1 = sum[i] - sum[st];
    //double sumsq1 = sumsq[i] - sumsq[st];

    double sum1 = sum[i];             
    double sumsq1 = sumsq[i];         
    if (buf_mid > w_length) {
        sum1 -= sum[st];    
        sumsq1 -= sumsq[st];
    }

    float sum2 = (float)(sum[en] - sum[i]);
    float sumsq2 = (float)(sumsq[en] - sumsq[i]);
    float mean1 = sum1 / w_lengthf;
    float mean2 = sum2 / w_lengthf;
    float combined_var = sumsq1 / w_lengthf - mean1 * mean1
        + sumsq2 / w_lengthf - mean2 * mean2;

    // Prevent problem due to very small variances
    combined_var = fmaxf(combined_var, eta);

    //t-stat
    //  Formula is a simplified version of Student's t-statistic for the
    //  special case where there are two samples of equal size with
    //  differing variance
    const float delta_mean = mean2 - mean1;
    return fabs(delta_mean) / sqrt(combined_var / w_lengthf);
}

bool EventDetector::peak_detect(float current_value, Detector &detector) {

    //Carry on if we've been masked out
    if (detector.masked_to >= buf_mid) {
        return false;
    }

    //float current_value = raw_buf[buf_center()];

    if (detector.peak_pos == detector.DEF_PEAK_POS) {
        //CASE 1: We've not yet recorded a maximum
        if (current_value < detector.peak_value) {
            //Either record a deeper minimum...
            detector.peak_value = current_value;
        } else if (current_value - detector.peak_value >
                   PRMS.peak_height) {
            // ...or we've seen a qualifying maximum
            detector.peak_value = current_value;
            detector.peak_pos = buf_mid;
            //otherwise, wait to rise high enough to be considered a peak
        }
    } else {
        //CASE 2: In an existing peak, waiting to see if it is good
        if (current_value > detector.peak_value) {
            //Update the peak
            detector.peak_value = current_value;
            detector.peak_pos = buf_mid;
        }
        //Dominate other tstat signals if we're going to fire at some point
        if (detector.window_length == short_detector.window_length) {
            if (detector.peak_value > detector.threshold) {
                long_detector.masked_to =
                    detector.peak_pos + detector.window_length;
                long_detector.peak_pos =
                    long_detector.DEF_PEAK_POS;
                long_detector.peak_value =
                    long_detector.DEF_PEAK_VAL;
                long_detector.valid_peak = false;
            }
        }
        //Have we convinced ourselves we've seen a peak
        if (detector.peak_value - current_value > PRMS.peak_height
            && detector.peak_value > detector.threshold) {
            detector.valid_peak = true;
        }
        //Finally, check the distance if this is a good peak
        if (detector.valid_peak
            && (buf_mid - detector.peak_pos) >
            detector.window_length / 2) {
            detector.peak_pos = detector.DEF_PEAK_POS;
            detector.peak_value = current_value;
            detector.valid_peak = false;

            return true;
        }
    }

    return false;
}



/**  Create an event given boundaries
 *
 *   Note: Bounds are CADLAG (i.e. lower bound is contained in the interval but
 *   the upper bound is not).
 *
 *  @param start Index of lower bound
 *  @param end Index of upper bound
 *  @param sums
 *  @param sumsqs
 *  @param nsample  Total number of samples in read
 *
 *  @returns An initialised event.  A 'null' event is returned on error.
 **/
Event EventDetector::create_event(u32 evt_en) {
    //Event event = { 0 };

    u32 evt_en_buf = evt_en % BUF_LEN;

    event_.start = evt_st;
    event_.length = (float)(evt_en - evt_st);
    event_.mean = (sum[evt_en_buf] - evt_st_sum) / event_.length;
    //const float deltasqr = (sumsq[evt_en_buf] - evt_st_sumsq);
    //const float var = deltasqr / event_.length - event_.mean * event_.mean;
    //event_.stdv = sqrtf(fmaxf(var, 0.0f));

    event_.mean = calibrate(event_.mean);
    //event_.stdv = calibrate(event_.stdv);

    evt_st = evt_en;
    evt_st_sum = sum[evt_en_buf];
    evt_st_sumsq = sumsq[evt_en_buf];

    len_sum_ += event_.length;
    total_events_++;

    return event_;
}


//EventDetectorMeta::EventDetectorMeta(Params prms) : EventDetector(prms) {
//    
//}

/*
EventDetectorMeta::EventDetectorMeta() : EventDetectorMeta(PRMS_DEF) {
    dbg_.init_window(long_detector.window_length);
    dbg_.init_window(short_detector.window_length);
}

void EventDetectorMeta::reset() {
    EventDetector::reset();
    dgb_.clear();
}

float EventDetectorMeta::compute_tstat(u32 w_length) {
    auto tstat = EventDetector::compute_tstat(w_length);
    dbg_[w_length].tstats.push_back(tstat);
    return tstat;
}

bool EventDetectorMeta::peak_detect(float current_value, EventDetector::Detector &detector) {
    auto peak = EventDetector::peak_detect(current_value, detector);
    if (peak) dbg_[detector.window_length].push_back(buf_mid);
    return peak;
}
*/
//=====================stop===================

