#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stddef.h>
#include <list>
#include <set>
#include <utility>
#include <unordered_map>
#include <sdsl/suffix_arrays.hpp>
#include "timer.hpp"
#include "leaf_aligner.hpp"
#include "fmi.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)
//#define DEBUG_PATHS
//#define DEBUG_RANGES

#define MAX_CHILDREN 5

unsigned char LeafAligner::Alignment::PROB_WIN_LEN = 0, 
              LeafAligner::Alignment::TYPE_WIN_LEN = 0;

//Source constructor
LeafAligner::Alignment::Alignment(Kmer kmer, double prob)
    : prob_sums_(new double[PROB_WIN_LEN]),
      event_types_(new EventType[TYPE_WIN_LEN]) {
    init_source(kmer, prob);
}

//Sibling constructor
LeafAligner::Alignment::Alignment(Alignment *a, Kmer kmer, 
                                  double prob, EventType type) 
    : prob_sums_(new double[PROB_WIN_LEN]),
      event_types_(new EventType[TYPE_WIN_LEN]) {
    init_from_sibling(a, kmer, prob, type);
}

LeafAligner::Alignment::~Alignment() {
    delete[] prob_sums_;
    delete[] event_types_;
}

void LeafAligner::Alignment::init_source(Kmer kmer, double prob) {
    length_ = 1;
    consec_stays_ = 0;
    prtl_ = 0;
    prhd_ = 1;
    prlen_ = 1;
    tytl_ = 0;
    tyhd_ = 0;
    tylen_ = 1;
    prev_kmer_ = kmer;

    //MATCH should always be index 0
    all_type_counts_[EventType::MATCH] = 1;
    win_type_counts_[EventType::MATCH] = 1;
    for (unsigned char t = 1; t < EventType::NUM_TYPES; t++) {
        all_type_counts_[t] = 0;
        win_type_counts_[t] = 0;
    }

    prob_sums_[prtl_] = 0;
    prob_sums_[prhd_] = prob;
    event_types_[tyhd_] = EventType::MATCH;
}

void LeafAligner::Alignment::init_from_sibling(Alignment *a, Kmer kmer, 
                                               double prob, EventType type) {
    length_ = a->length_;
    prtl_ = a->prtl_;
    prhd_ = a->prhd_;
    prlen_ = a->prlen_;
    tytl_ = a->tytl_;
    tyhd_ = a->tyhd_;
    tylen_ = a->tylen_;
    prev_kmer_ = kmer;

    std::memcpy(prob_sums_, a->prob_sums_, PROB_WIN_LEN * sizeof(double));
    std::memcpy(event_types_, a->event_types_, TYPE_WIN_LEN * sizeof(EventType));
    std::memcpy(all_type_counts_, a->all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, a->win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));


    prob_sums_[prhd_] = 
        prob_sums_[prhd_ != 0 ? prhd_ - 1 : PROB_WIN_LEN - 1] + prob;

    all_type_counts_[event_types_[tyhd_]]--;
    all_type_counts_[type]++;

    win_type_counts_[event_types_[tyhd_]]--;
    win_type_counts_[type]++;

    event_types_[tyhd_] = type;

    //"Sibling" should still have parent's conesc_stays_
    consec_stays_ = type == EventType::STAY ? a->consec_stays_ + 1 : 0;
}

void LeafAligner::Alignment::init_from_parent(Alignment *a, Kmer kmer, 
                                              double prob, EventType type) {
    length_ = a->length_;
    prtl_ = a->prtl_;
    prhd_ = a->prhd_;
    prlen_ = a->prlen_;
    tytl_ = a->tytl_;
    tyhd_ = a->tyhd_;
    tylen_ = a->tylen_;
    consec_stays_ = a->consec_stays_;
    prev_kmer_ = kmer;

    std::memcpy(prob_sums_, a->prob_sums_, PROB_WIN_LEN * sizeof(double));
    std::memcpy(event_types_, a->event_types_, TYPE_WIN_LEN * sizeof(EventType));
    std::memcpy(all_type_counts_, a->all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, a->win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));

    make_child(kmer, prob, type);
}

void LeafAligner::Alignment::make_child(Kmer kmer, 
                                        double prob, 
                                        EventType type) {
    length_++;
    prev_kmer_ = kmer;

    if (tylen_ < TYPE_WIN_LEN) { //TODO: maybe off by 1
        tyhd_++;
        tylen_++;
    } else {
        //std::cout << "incr type " << (int) tyhd_ << " " << event_types_[tyhd_] << " " << (int) win_type_counts_[type] << " ";

        tyhd_ = tytl_;
        tytl_ = (tytl_ + 1) % TYPE_WIN_LEN;
        win_type_counts_[event_types_[tyhd_]]--;

        win_type_counts_[event_types_[tytl_]]--;
        event_types_[tytl_] = EventType::MATCH;
        win_type_counts_[EventType::MATCH]++;

        //std::cout << (int) tyhd_ << " " << event_types_[tyhd_]  << " " << (int) win_type_counts_[type] << " ";
    }
    //std::cout << (int) win_type_counts_[type] << "\n";

    event_types_[tyhd_] = type;
    win_type_counts_[type]++;
    all_type_counts_[type]++;

    if (prlen_ < PROB_WIN_LEN - 1) { //TODO: maybe off by 1
        prob_sums_[prhd_+1] = prob_sums_[prhd_] + prob;
        prhd_++;
        prlen_++;
    } else {
        prob_sums_[prtl_] = prob_sums_[prhd_] + prob;
        prhd_ = prtl_;
        prtl_ = (prtl_ + 1) % PROB_WIN_LEN;
        //prlen_ = PROB_WIN_LEN;
    }

    //Still need to update consec stays
}

//Should be called on first child after all siblings have been created
void LeafAligner::Alignment::update_consec_stays() {
    if (event_types_[tyhd_] == EventType::STAY) {
        consec_stays_++;
    } else {
        consec_stays_ = 0;
    }
}

size_t LeafAligner::Alignment::event_len() {
    return length_;
}

size_t LeafAligner::Alignment::match_len() {
    return win_type_counts_[EventType::MATCH];
}

double LeafAligner::Alignment::mean_prob() const {
    return (prob_sums_[prhd_] - prob_sums_[prtl_]) / prlen_;
}

//double LeafAligner::Alignment::next_mean_prob(double next_prob) const {
//    return (seed_prob_ + next_prob) / (length_ + 1);
//}

bool LeafAligner::Alignment::better_than_sibling(const Alignment *a, double prob) {
    double x = a->prob_sums_[a->prhd_ > 0 ? a->prhd_ - 1 : PROB_WIN_LEN-1],
           y = a->prob_sums_[a->prtl_];

    double replace_prob = (x - y + prob) / a->prlen_;
    //std::cout << (int) a->prtl_ << "\t" << (int) a->prhd_ << "\t"
    //          << a->prob_sums_[a->prtl_] << "\t"
    //          << a->prob_sums_[a->prhd_] << "\t" << x << "\t"
    //          << replace_prob << "\t" << mean_prob() << "\n";

    return replace_prob > mean_prob();
}

bool LeafAligner::Alignment::better_than_parent(const Alignment *a, double prob) {
    double x = a->prob_sums_[a->prhd_], 
           y = a->prob_sums_[(a->prlen_ == PROB_WIN_LEN-1 && a->prtl_ < PROB_WIN_LEN-1) ? a->prtl_ + 1 : 0];
    unsigned int len = a->prlen_ + (a->prlen_ < PROB_WIN_LEN-1);
    double replace_prob = (x - y + prob) / len;
    //std::cout << (int) a->prtl_ << "\t" << (int) a->prhd_ << "\n";

    return replace_prob > mean_prob();
}

double LeafAligner::Alignment::next_mean_prob() {
    double x = prob_sums_[prhd_], 
           y = prob_sums_[(prlen_ == PROB_WIN_LEN-1 && prtl_ < PROB_WIN_LEN-1) ? prtl_ + 1 : 0];
    unsigned int len = prlen_ - (prlen_ == PROB_WIN_LEN-1);
    return (x - y) / len;

}

bool LeafAligner::Alignment::should_report(const AlnParams &params) {
    bool r = event_types_[tyhd_] == EventType::MATCH && 
           win_type_counts_[EventType::STAY] <= params.max_stay_frac_ * tylen_ &&
           mean_prob() >= params.min_seed_pr_;
    //std::cout << (int) win_type_counts_[EventType::STAY] << " stays " << (int) tylen_ << " " << r << "\n";
    return r;
}

LeafAligner::LeafAligner(const FMI &fmi, 
                     const AlnParams &ap,
                     const std::string &label)
    : fmi_(fmi),
      params_(ap),
      label_(label) {
    timer.reset();

    Alignment::PROB_WIN_LEN = params_.graph_elen_+1;
    Alignment::TYPE_WIN_LEN = params_.graph_elen_;

    kmer_ranges_ = new Range[params_.model_.kmer_count()];
    for (Kmer k = 0; k < params_.model_.kmer_count(); k++) {
        Range r = fmi_.get_full_range(params_.model_.get_last_base(k));
        for (size_t i = params_.model_.kmer_len()-2; 
             i < params_.model_.kmer_len(); i--) {
            r = fmi_.get_neighbor(r, params_.model_.get_base(k, i));
        }
        kmer_ranges_[k] = r;
    }
}

LeafAligner::~LeafAligner() {
    delete[] kmer_ranges_;
    reset();
}

void LeafAligner::new_read(size_t read_len) {
    reset();
    cur_event_ = read_len;
}

void LeafAligner::reset() {
    for (auto p = prev_alns_.begin(); p != prev_alns_.end(); p++) {
        inactive_alns_.push_back(p->second);
    }
    prev_alns_.clear();
}

std::vector<Result> LeafAligner::add_event(double *kmer_probs, std::ostream &out) {

    //Update event index
    cur_event_--;

    Timer t;

    double prob;

    t.reset();

    Range prev_range;
    Alignment *prev_aln;
    Kmer prev_kmer;
    
    bool child_found;
    double evpr_thresh;


    //Find neighbors of previous nodes
    for (auto p = prev_alns_.begin(); p != prev_alns_.end(); p++) {
        Range prev_range = p->first;
        prev_aln = p->second;


        //std::cout << cur_event_ << "\t"
        //          << prev_range.start_ << "-" << prev_range.end_ << "\t"
        //          << (int) prev_aln->win_type_counts_[EventType::STAY] << "\t"
        //          << prev_aln->length_ << "\t"
        //          << (int) prev_aln->prlen_ << "\t"
        //          << prev_aln->next_mean_prob() << "\n";


        child_found = false;

        evpr_thresh = 0;
        for (size_t i = 0; i < params_.expr_lengths_.size(); i++) {
            if (prev_range.length() <= params_.expr_lengths_[i]) {
                evpr_thresh = params_.expr_probs_[i];
                break;
            }
        }

        if (evpr_thresh == 0) {
            evpr_thresh = params_.min_anchor_evpr_;
        }

        //Get probability for stay neighbor
        prev_kmer = prev_aln->prev_kmer_;
        prob = kmer_probs[prev_kmer];

        if (prev_aln->consec_stays_ < params_.max_consec_stay_ && prob >= evpr_thresh) {
            //std::cout << prev_range.start_ << "-" << prev_range.end_ << "\t" << prob << "\n";
            child_found = 
                add_child(prev_range,
                          prev_aln, 
                          prev_kmer, 
                          prob, 
                          EventType::STAY, 
                          child_found);
        }
        
        //Find next possible kmers
        auto neighbor_itr = params_.model_.get_neighbors(prev_kmer);
        std::list<Kmer> next_kmers;

        for (auto n = neighbor_itr.first; n != neighbor_itr.second; n++) {
            prob = kmer_probs[*n];
            if(prob >= evpr_thresh) { //&& prev_node->next_mean_prob(prob) >= params_.min_seed_pr_) {
                next_kmers.push_back(*n);
            } 
        }

        //Find ranges FM index for those next kmers
        //std::list<Range> next_ranges = fmi_.get_neigbhors(prev_range, next_kmers);

        auto next_kmer = next_kmers.begin(); 
        //auto next_range = next_ranges.begin();
        

        //Add all the neighbors that were found
        while (next_kmer != next_kmers.end()) {
            prob = kmer_probs[*next_kmer];
            //rank = kmer_ranks[*next_kmer];

            Base next_base = params_.model_.get_first_base(*next_kmer);
            Range next_range = fmi_.get_neighbor(prev_range, next_base);

            //std::cout << next_range.start_ << "-" << next_range.end_ << "\t" << prob << "\n";

            child_found = 
                add_child(next_range,
                          prev_aln, 
                          *next_kmer, 
                          prob, 
                          EventType::MATCH, 
                          child_found);


            next_kmer++; 
            //next_range++;
        }

        if (child_found) {
            prev_aln->update_consec_stays();
        } else {
            inactive_alns_.push_back(prev_aln);
        }
    }

    //Find sources
    for (Kmer kmer = 0; kmer < params_.model_.kmer_count(); kmer++) {
        prob = kmer_probs[kmer];
        if (prob >= params_.min_anchor_evpr_) {
            Range next_range = kmer_ranges_[kmer];

            if (next_range.is_valid()) {
                //Will split sources if they intersect existing nodes
                add_sources(next_range, kmer, prob);             
            }
        }
    }

    prev_alns_.clear();
    prev_alns_.swap(next_alns_);

    auto r = pop_seeds(out);

    return r;
}

size_t LeafAligner::add_sources(const Range &range, Kmer kmer, double prob) {

    //Find closest existing node
    auto lb = next_alns_.lower_bound(range);
    
    //Find range of existing alns intersecting the new aln
    auto start = lb;
    while (start != next_alns_.begin() 
           && std::prev(start)->first.intersects(range)){
        start--;
    }

    auto end = lb;
    while (end != next_alns_.end() && end->first.intersects(range)) {
        end++;
    }

    //No alns intersect new aln, just add it
    if (start == next_alns_.end() || (start == lb && end == lb)) {

        if (inactive_alns_.empty()) {
            next_alns_[range] = new Alignment(kmer, prob);
        } else {
            next_alns_[range] = inactive_alns_.back();
            next_alns_[range]->init_source(kmer, prob);
            inactive_alns_.pop_back();
        }

        return 1;
    }

    //Split range to not include intersecting ranges
    std::vector<Range> split_ranges;
    //TODO: pre-alocate based on start/end distance?

    Range rr(range); //Copy full range

    for (auto pr = start; pr != end; pr++) {

        //After this, rr stores right segment, lr stores left
        Range lr = rr.split_range(pr->first); //Left range

        //Add left range if it exists
        if (lr.is_valid()) {
            split_ranges.push_back(lr);
        }
        
        //If right range invalid, no more splitting can be done
        if (!rr.is_valid()) {
            break;
        }
    }

    //Add last range if there is one
    if (rr.is_valid()) {
        split_ranges.push_back(rr);
    }

    //Add a new aln for every split range
    for (auto r = split_ranges.begin(); r != split_ranges.end(); r++) {
        if (inactive_alns_.empty()) {
            next_alns_[*r] = new Alignment(kmer, prob);
        } else {
            next_alns_[*r] = inactive_alns_.back();
            next_alns_[*r]->init_source(kmer, prob);
            inactive_alns_.pop_back();
        }
    }

    return split_ranges.size();
}

bool LeafAligner::add_child(Range &range, 
                            Alignment *prev_aln,
                            Kmer kmer,
                            double prob,
                            EventType type,
                            bool prev_is_sibling) {
    
    if (!range.is_valid()) {
        return prev_is_sibling;
    }
    
    //Find closest node >= the node being added
    auto lb = next_alns_.lower_bound(range);

    //std::cout << "a\n";

    //Alignment range hasn't been added yet
    if (lb == next_alns_.end() || !lb->first.same_range(range)) {

        Alignment *new_aln;
        if (prev_is_sibling) {
            if (inactive_alns_.empty()) {
                //std::cout << "b\n";
                new_aln = new Alignment(prev_aln, kmer, prob, type);
            } else {
                //std::cout << "c\n";
                new_aln = inactive_alns_.back();
                inactive_alns_.pop_back();
                new_aln->init_from_sibling(prev_aln, kmer, prob, type);
            }
        } else {
            //std::cout << "c/b\n";
            prev_aln->make_child(kmer, prob, type);
            new_aln = prev_aln;
        }

        //Store it with it's range
        next_alns_.insert(lb, std::pair<Range, Alignment *>(range, new_aln));

        //Created one node
        return true;
    }
    //std::cout << "old... \n";

    
    //Alignment associated with same range
    Alignment *dup_aln = lb->second;

              
    //std::cout << "d " << prev_is_sibling << "\n";

    if (prev_is_sibling) {
        if (dup_aln->better_than_sibling(prev_aln, prob)) {
            //std::cout << "e1\n";
            dup_aln->init_from_sibling(prev_aln, kmer, prob, type);
            return true;
        }
    } else {
        if (dup_aln->better_than_parent(prev_aln, prob)) {
            //std::cout << "e2\n";
            inactive_alns_.push_back(lb->second);
            lb->second = prev_aln;
            lb->second->make_child(kmer, prob, type);
            return true;
        }
    }
    //std::cout << "f\n";

    return prev_is_sibling;
}



std::vector<Result> LeafAligner::pop_seeds(std::ostream &out) { 

    std::vector<Result> results;

    for (auto p = prev_alns_.begin(); p != prev_alns_.end(); p++) {

        Alignment *n = p->second;
        //std::cout << n->event_len() << "\t" << params_.graph_elen_ << "\n";

        if (n->event_len() >= params_.graph_elen_) {

            Range range = p->first;

            //TODO: max repeat parameter
            if (n->should_report(params_) && range.length() < 100) { 
                Result r(cur_event_, params_.graph_elen_, n->mean_prob());

                for (unsigned int s = range.start_; s <= range.end_; s++) {
                    r.set_ref_range(fmi_.sa(s), n->match_len());
                    results.push_back(r);

                    out << label_ << "\t";

                    r.print(out);
                }
            }
        }
    }
        

    return results;

}

