#ifndef ALIGNMENT_STRUCTURE_HPP
#define ALIGNMENT_STRUCTURE_HPP

#include "fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"
#include <list>
#include <iostream>

//#define VERBOSE_TIME

class Result {
    public:

    Result(unsigned int read_start, 
           unsigned int seed_len, 
           float prob, 
           unsigned int ref_start = 0, 
           unsigned int ref_end = 0);

    void set_ref_range(unsigned int start, unsigned int end);
    void print(std::ostream &out);

    Range read_range_, ref_range_;
    float seed_prob_;

    #ifdef DEBUG_PROB
    float min_evt_prob_;
    #endif
};

class AlnParams {
    public:
    const KmerModel &model_;

    unsigned int path_win_len_, 
                 min_rep_len_,
                 max_rep_copy_,
                 max_paths_;

    float max_stay_frac_;

    unsigned int max_consec_stay_, 
                 max_ignores_,
                 max_skips_;

    float window_prob_;

    std::vector<unsigned int> evpr_lengths_;
    std::vector<float> evpr_threshes_;

    AlnParams(const KmerModel &model,
              unsigned int path_win_len, 
              unsigned int min_rep_len, 
              unsigned int max_rep_copy, 
              unsigned int max_paths, 
              float max_stay_frac,
              unsigned int max_consec_stay,
              unsigned int max_ignores, 
              unsigned int max_skips,
              const std::string event_probs,
              float window_prob);           
    
    unsigned int nucl_to_events(unsigned int n);
    float get_prob_thresh(unsigned int fm_length);
    float get_source_prob();
};


class Aligner {

    public:

    virtual void new_read(size_t read_len) = 0;
    virtual void reset() = 0;
    virtual std::vector<Result> add_event(float *kmer_probs, std::ostream &seeds_out, std::ostream &time_out) = 0;
};

#endif
