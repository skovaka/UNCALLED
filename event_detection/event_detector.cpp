#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>

#include "event_detector.hpp"
extern "C" {
    #include "scrappie_stdlib.h"
}

typedef Detector *DetectorPtr;

//window_length1
//window_length2
//threshold1
//threshold2
//peak_height

EventDetector::EventDetector(const detector_param &edparams) :
    params (edparams),
    BUF_LEN (1 + edparams.window_length2 * 2) {

    //event_table et = { 0 };

    //sums = (double *) calloc(rt.n + 1, sizeof(double));
    //sumsqs = (double *) calloc(rt.n + 1, sizeof(double));

    sum = new double[BUF_LEN];
    sumsq = new double[BUF_LEN];

    sum[0] = sumsq[0] = 0.0;
    t = t_mod = 1;
    evt_st = 0;
    evt_st_sum = evt_st_sumsq = 0.0;

    short_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .threshold = edparams.threshold1,
        .window_length = edparams.window_length1,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false
    };

    
                               


    long_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .threshold = edparams.threshold2,
        .window_length = edparams.window_length2,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false
    };
}

size_t EventDetector::sig_center() {
    return t - (BUF_LEN / 2) - 1;
}

event_t EventDetector::add_signal(float s) {
    
    if (t_mod > 0) {
        sum[t_mod] = sum[t_mod-1] + s;
        sumsq[t_mod] = sumsq[t_mod-1] + s*s;
    } else {
        sum[t_mod] = sum[BUF_LEN-1] + s;
        sumsq[t_mod] = sumsq[BUF_LEN-1] + s*s;
    }


    t++;
    if (++t_mod == BUF_LEN) {
        t_mod = 0;
    }


    event_t ret = {0};

    //if (t >= BUF_LEN) {
    int sc = sig_center();
    double tstat1 = compute_tstat(sc, params.window_length1),
           tstat2 = compute_tstat(sc, params.window_length2);


    bool p1 = peak_detect(tstat1, short_detector),
         p2 = peak_detect(tstat2, long_detector);

    if (p1 || p2) {
        ret = create_event();
    }
    //}
    
    return ret;
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
float EventDetector::compute_tstat(size_t s, size_t w_length) {
    assert(w_length > 0);

    //float *tstat = (float *) calloc(d_length, sizeof(float));

    const float eta = FLT_MIN;
    const float w_lengthf = (float) w_length;

    // Quick return:
    //   t-test not defined for number of points less than 2
    //   need at least as many points as twice the window length
    if (t <= 2 * w_length || w_length < 2) {
        return 0;
    }

    // fudge boundaries
    //for (size_t i = 0; i < w_length; ++i) {
    //    tstat[i] = 0;
    //    tstat[d_length - i - 1] = 0;
    //}

    size_t i = s % BUF_LEN,
           st = (s - w_length) % BUF_LEN,
           en = (s + w_length) % BUF_LEN;

    //std::cout << i << " " << st << " " << en << "\n";

    double sum1 = sum[i];
    double sumsq1 = sumsq[i];
    //if (i > w_length) {
        sum1 -= sum[st];
        sumsq1 -= sumsq[st];
    //}
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

    size_t i = sig_center();

    printf("s %d\t%d\t%f\t%d\n", 
           short_detector.masked_to,
           short_detector.peak_pos,
           short_detector.peak_value,
           short_detector.valid_peak);
    printf("l %d\t%d\t%f\t%d\n", 
           long_detector.masked_to,
           long_detector.peak_pos,
           long_detector.peak_value,
           long_detector.valid_peak);

    //Carry on if we've been masked out
    if (detector.masked_to >= i) {
        return false;
    }

    //float current_value = raw_buf[buf_center()];

    if (detector.peak_pos == detector.DEF_PEAK_POS) {
        //CASE 1: We've not yet recorded a maximum
        if (current_value < detector.peak_value) {
            //Either record a deeper minimum...
            detector.peak_value = current_value;
        } else if (current_value - detector.peak_value >
                   params.peak_height) {
            // ...or we've seen a qualifying maximum
            detector.peak_value = current_value;
            detector.peak_pos = i;
            //otherwise, wait to rise high enough to be considered a peak
        }
    } else {
        //CASE 2: In an existing peak, waiting to see if it is good
        if (current_value > detector.peak_value) {
            //Update the peak
            detector.peak_value = current_value;
            detector.peak_pos = i;
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
        if (detector.peak_value - current_value > params.peak_height
            && detector.peak_value > detector.threshold) {
            detector.valid_peak = true;
        }
        //Finally, check the distance if this is a good peak
        if (detector.valid_peak
            && (i - detector.peak_pos) >
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
event_t EventDetector::create_event() {
    event_t event = { 0 };
    event.pos = -1;
    event.state = -1;

    int evt_en = sig_center(),
        evt_en_buf = evt_en % BUF_LEN;


    event.start = (uint64_t) evt_st;
    event.length = (float)(evt_en - evt_st); //TODO: OFF BY ONE?
    event.mean = (float)(sum[evt_en_buf] - evt_st_sum) / event.length;
    const float deltasqr = (sumsq[evt_en_buf] - evt_st_sumsq);
    const float var = deltasqr / event.length - event.mean * event.mean;
    event.stdv = sqrtf(fmaxf(var, 0.0f));

    evt_st = evt_en;
    evt_st_sum = sum[evt_en_buf];
    evt_st_sumsq = sumsq[evt_en_buf];

    return event;
}

//=====================stop===================

