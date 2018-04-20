#include "aligner.hpp"

AlnParams::AlnParams(const KmerModel &model,
                     unsigned int path_win_len, 
                     unsigned int min_rep_len, 
                     unsigned int max_rep_copy, 
                     double max_stay_frac,
                     unsigned int max_consec_stay,
                     unsigned int max_ignores, 
                     unsigned int max_skips,
                     const std::string event_probs,
                     double window_prob) 
        : model_(model),
          path_win_len_(path_win_len),
          min_rep_len_(min_rep_len),
          max_rep_copy_(max_rep_copy),
          max_stay_frac_(max_stay_frac),
          max_consec_stay_(max_consec_stay),
          max_ignores_(max_ignores),
          max_skips_(max_skips),
          window_prob_(window_prob) {
    
    //"-2.25_100-2.4_5-4.0_1-10.0"
    //1-10.0_5-4.0_100-2.4
    size_t i = event_probs.find('_'), j, k;

    evpr_threshes_.push_back(atof(event_probs.substr(0, i).c_str()));

    i += 1;
    j = event_probs.find('_');

    while(i < event_probs.size()) {
        k = event_probs.find('-', i);
        evpr_lengths_.push_back( atoi(event_probs.substr(i, k).c_str()) );
        evpr_threshes_.push_back( atof(event_probs.substr(k, j).c_str()) ); 

        i = j+1;
        j = event_probs.find('_', i+1);
        if (j == std::string::npos) {
            j = event_probs.size();
        }
    }
}

double AlnParams::get_prob_thresh(unsigned int fm_length) {
    auto pr = evpr_threshes_.begin();
    for (auto len = evpr_lengths_.begin(); len != evpr_lengths_.end(); len++) {
        if (fm_length > *len) {
            break;
        }
        pr++;
    }

    //std::cout << fm_length << "\t" << (*pr) << "\n";

    return *pr;
}

double AlnParams::get_source_prob() {
    return evpr_threshes_.front();
}

unsigned int AlnParams::nucl_to_events(unsigned int n) {
    return n - model_.kmer_len() + 1;
}

//unsigned int AlnParams::get_graph_len(unsigned int seed_nlen) {
//    return (nucl_to_events(seed_nlen) / (1.0 - max_stay_frac_)) + max_ignores_;
//}

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
