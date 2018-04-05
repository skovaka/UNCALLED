#include "aligner.hpp"

AlnParams::AlnParams(const KmerModel &model,
                     unsigned int min_seed_nlen, 
                     unsigned int anchor_nlen, 
                     unsigned int max_ignores, 
                     unsigned int max_skips,
                      unsigned int max_consec_stay,
                     double max_stay_frac,
                     double min_anchor_evpr,
                     //double min_extend_evpr,
                     std::vector<unsigned int>    expr_lengths,
                     std::vector<double> expr_probs,
                     double min_seed_pr,
                     double min_stay_pr)
        : model_(model),
          max_ignores_(max_ignores),
          max_skips_(max_skips),
          max_consec_stay_(max_consec_stay),
          max_stay_frac_(max_stay_frac),
          min_anchor_evpr_(min_anchor_evpr),
          min_seed_pr_(min_seed_pr),
          min_stay_pr_(min_stay_pr),
          //min_extend_evpr_(min_extend_evpr),
          expr_lengths_(expr_lengths),
          expr_probs_(expr_probs) {

    anchor_rlen_ = nucl_to_events(anchor_nlen);
    graph_elen_ = get_graph_len(min_seed_nlen);

    std::cerr << "Graph len: " << graph_elen_ << "\n";

}

unsigned int AlnParams::nucl_to_events(unsigned int n) {
    return n - model_.kmer_len() + 1;
}

unsigned int AlnParams::get_graph_len(unsigned int seed_nlen) {
    return (nucl_to_events(seed_nlen) / (1.0 - max_stay_frac_)) + max_ignores_;
}

Result::Result(unsigned int read_start, 
              unsigned  int seed_len, 
              double prob, 
              unsigned int ref_start, 
              unsigned int ref_end) 
    : read_range_( Range(read_start, read_start + seed_len - 1) ),
      ref_range_(Range(ref_start, ref_end)),
      seed_prob_(prob) {}

void Result::set_ref_range(unsigned int start, unsigned int length) {
    ref_range_.start_ = start;
    ref_range_.end_ = start + length - 1;
}


void Result::print(std::ostream &out) {
    out << read_range_.start_ << "-" << read_range_.end_ << "\t"
              << ref_range_.start_ << "-" << ref_range_.end_ << "\t"
              << seed_prob_;
              

    #ifdef DEBUG_PROB
    out << " " << min_evt_prob_;
    #endif

    out << "\n";
}
