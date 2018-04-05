#ifndef ALIGNMENT_STRUCTURE_HPP
#define ALIGNMENT_STRUCTURE_HPP

#include "fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"
#include <list>
#include <iostream>

class Result {
    public:

    Result(unsigned int read_start, 
           unsigned int seed_len, 
           double prob, 
           unsigned int ref_start = 0, 
           unsigned int ref_end = 0);

    void set_ref_range(unsigned int start, unsigned int end);
    void print(std::ostream &out);

    Range read_range_, ref_range_;
    double seed_prob_;

    #ifdef DEBUG_PROB
    double min_evt_prob_;
    #endif
};

class AlnParams {
    public:
    const KmerModel &model_;

    size_t graph_elen_, anchor_rlen_, max_ignores_, max_skips_, max_consec_stay_;

    double max_stay_frac_,
           min_anchor_evpr_,
           min_seed_pr_,
           min_stay_pr_;
           //min_extend_evpr_;

    std::vector<unsigned int> expr_lengths_;
    std::vector<double> expr_probs_;

    AlnParams(const KmerModel &model,
              unsigned int min_seed_nlen, 
              unsigned int anchor_nlen, 
              unsigned int max_ignores, 
              unsigned int max_skips,
              unsigned int max_consec_stay,
              double max_stay_frac,
              double min_anchor_evpr,
              //double min_extend_evpr,
              std::vector<unsigned int> expr_lengths,
              std::vector<double> expr_probs,
              double min_seed_pr,
              double min_stay_pr);
    
    unsigned int nucl_to_events(unsigned int n);
    unsigned int get_graph_len(unsigned int seed_nlen);
};


class Aligner {

    public:

    virtual void new_read(size_t read_len) = 0;
    virtual void reset() = 0;
    virtual std::vector<Result> add_event(Event e, std::ostream &out) = 0;
};

#endif
