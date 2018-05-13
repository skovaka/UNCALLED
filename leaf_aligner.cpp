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

unsigned char LeafAligner::PathBuffer::PROB_WIN_LEN = 0, 
              LeafAligner::PathBuffer::TYPE_WIN_LEN = 0;

//Source constructor
LeafAligner::PathBuffer::PathBuffer(Kmer kmer, float prob)
    : prob_sums_(new float[PROB_WIN_LEN]),
      event_types_(new EventType[TYPE_WIN_LEN]) {
    init_source(kmer, prob);
}

//Sibling constructor
LeafAligner::PathBuffer::PathBuffer(PathBuffer *a, Kmer kmer, 
                                  float prob, EventType type) 
    : prob_sums_(new float[PROB_WIN_LEN]),
      event_types_(new EventType[TYPE_WIN_LEN]) {
    init_from_sibling(a, kmer, prob, type);
}

LeafAligner::PathBuffer::~PathBuffer() {
    delete[] prob_sums_;
    delete[] event_types_;
}

void LeafAligner::PathBuffer::init_source(Kmer kmer, float prob) {
    length_ = 1;
    consec_stays_ = 0;
    prtl_ = 0;
    prhd_ = 1;
    prlen_ = 1;
    tytl_ = 0;
    tyhd_ = 0;
    tylen_ = 1;
    prev_kmer_ = kmer;
    sa_checked_ = false;

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

void LeafAligner::PathBuffer::init_from_sibling(PathBuffer *a, Kmer kmer, 
                                               float prob, EventType type) {
    length_ = a->length_;
    prtl_ = a->prtl_;
    prhd_ = a->prhd_;
    prlen_ = a->prlen_;
    tytl_ = a->tytl_;
    tyhd_ = a->tyhd_;
    tylen_ = a->tylen_;
    prev_kmer_ = kmer;
    sa_checked_ = a->sa_checked_;

    std::memcpy(prob_sums_, a->prob_sums_, PROB_WIN_LEN * sizeof(float));
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

void LeafAligner::PathBuffer::init_from_parent(PathBuffer *a, Kmer kmer, 
                                              float prob, EventType type) {
    length_ = a->length_;
    prtl_ = a->prtl_;
    prhd_ = a->prhd_;
    prlen_ = a->prlen_;
    tytl_ = a->tytl_;
    tyhd_ = a->tyhd_;
    tylen_ = a->tylen_;
    consec_stays_ = a->consec_stays_;
    prev_kmer_ = kmer;
    sa_checked_ = a->sa_checked_;


    std::memcpy(prob_sums_, a->prob_sums_, PROB_WIN_LEN * sizeof(float));
    std::memcpy(event_types_, a->event_types_, TYPE_WIN_LEN * sizeof(EventType));
    std::memcpy(all_type_counts_, a->all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, a->win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));

    make_child(kmer, prob, type);
}

void LeafAligner::PathBuffer::make_child(Kmer kmer, 
                                        float prob, 
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
void LeafAligner::PathBuffer::update_consec_stays() {
    if (event_types_[tyhd_] == EventType::STAY) {
        consec_stays_++;
    } else {
        consec_stays_ = 0;
    }
}

size_t LeafAligner::PathBuffer::event_len() {
    return length_;
}

size_t LeafAligner::PathBuffer::match_len() {
    return win_type_counts_[EventType::MATCH];
}

float LeafAligner::PathBuffer::mean_prob() const {
    return (prob_sums_[prhd_] - prob_sums_[prtl_]) / prlen_;

}

//float LeafAligner::PathBuffer::next_mean_prob(float next_prob) const {
//    return (seed_prob_ + next_prob) / (length_ + 1);
//}

bool LeafAligner::PathBuffer::better_than_sibling(const PathBuffer *a, float prob) {
    float x = a->prob_sums_[a->prhd_ > 0 ? a->prhd_ - 1 : PROB_WIN_LEN-1],
           y = a->prob_sums_[a->prtl_];

    float replace_prob = (x - y + prob) / a->prlen_;
    //std::cout << (int) a->prtl_ << "\t" << (int) a->prhd_ << "\t"
    //          << a->prob_sums_[a->prtl_] << "\t"
    //          << a->prob_sums_[a->prhd_] << "\t" << x << "\t"
    //          << replace_prob << "\t" << mean_prob() << "\n";

    return replace_prob > mean_prob();
}

bool LeafAligner::PathBuffer::better_than_parent(const PathBuffer *a, float prob) {
    float x = a->prob_sums_[a->prhd_], 
           y = a->prob_sums_[(a->prlen_ == PROB_WIN_LEN-1 && a->prtl_ < PROB_WIN_LEN-1) ? a->prtl_ + 1 : 0];
    unsigned int len = a->prlen_ + (a->prlen_ < PROB_WIN_LEN-1);
    float replace_prob = (x - y + prob) / len;
    //std::cout << (int) a->prtl_ << "\t" << (int) a->prhd_ << "\n";

    return replace_prob > mean_prob();
}

float LeafAligner::PathBuffer::next_mean_prob() {
    float x = prob_sums_[prhd_], 
           y = prob_sums_[(prlen_ == PROB_WIN_LEN-1 && prtl_ < PROB_WIN_LEN-1) ? prtl_ + 1 : 0];
    unsigned int len = prlen_ - (prlen_ == PROB_WIN_LEN-1);
    return (x - y) / len;

}

bool LeafAligner::PathBuffer::should_report(const Range &r, 
                                            const AlnParams &p,
                                            bool has_children) {
    return (r.length() == 1 || 
                (!has_children &&
                 r.length() <= p.max_rep_copy_ &&
                 match_len() >= p.min_rep_len_)) &&

           length_ >= p.path_win_len_ &&
           event_types_[tyhd_] == EventType::MATCH &&
           (!has_children || win_type_counts_[EventType::STAY] <= p.max_stay_frac_ * tylen_) &&
           mean_prob() >= p.window_prob_;
}

LeafAligner::LeafAligner(const FMI &fmi, 
                     const AlnParams &ap,
                     const std::string &label)
    : fmi_(fmi),
      params_(ap),
      label_(label) {

    PathBuffer::PROB_WIN_LEN = params_.path_win_len_+1;
    PathBuffer::TYPE_WIN_LEN = params_.path_win_len_;

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
    for (auto p = prev_paths_.begin(); p != prev_paths_.end(); p++) {
        inactive_paths_.push_back(p->second);
    }
    prev_paths_.clear();
}

std::vector<Result> LeafAligner::add_event(float *kmer_probs, 
                                           std::ostream &seeds_out,
                                           std::ostream &time_out) {

    //Update event index
    cur_event_--;

    float prob;

    Range prev_range;
    PathBuffer *prev_path;
    Kmer prev_kmer;
    
    bool child_found;
    float evpr_thresh;

    std::vector<Result> results;

    //std::cout << cur_event_ << "\t" << label_ << "\t" << prev_paths_.size() << "\n";
    //std::cout.flush()

    if (prev_paths_.size() > params_.max_paths_) {
        for (auto p = prev_paths_.begin(); p != prev_paths_.end(); p++) {
            inactive_paths_.push_back(p->second);
        }
        prev_paths_.clear();
    }

    #ifdef VERBOSE_TIME
    float stay_time = 0, 
           fm_time = 0, 
           neighbor_time = 0, 
           un_ext_time = 0;
    child_map_time_ = 0;
    child_add_time_ = 0;
    child_rpl_time_ = 0;
    Timer timer;
    #endif

    //Find neighbors of previous nodes
    for (auto p = prev_paths_.begin(); p != prev_paths_.end(); p++) {

        #ifdef VERBOSE_TIME
        //child_map_time_ += timer.get();
        //stay_time += timer.lap();
        #endif

        Range prev_range = p->first;
        prev_path = p->second;


        child_found = false;

        evpr_thresh = params_.get_prob_thresh(prev_range.length());

        //Get probability for stay neighbor
        prev_kmer = prev_path->prev_kmer_;
        prob = kmer_probs[prev_kmer];

        //std::cout << cur_event_ << "\t" 
        //          << prev_path->length_ << "\t"
        //          << prev_path->mean_prob() << "\t"
        //          << prev_kmer << "\t"
        //          << prev_path->event_types_[prev_path->tyhd_] << "\t"
        //          << prev_range.start_ << "-"
        //          << prev_range.end_ << "\n";

        #ifdef VERBOSE_TIME
        timer.reset();
        #endif

        if (prev_path->consec_stays_ < params_.max_consec_stay_ && 
            //prev_path->win_type_counts_[EventType::STAY] < params_.max_stay_frac_ * prev_path->tylen_ &&
            prob >= evpr_thresh) {
            child_found = 
                add_child(prev_range,
                          prev_path, 
                          prev_kmer, 
                          prob, 
                          EventType::STAY, 
                          child_found);
        }

        #ifdef VERBOSE_TIME
        stay_time += timer.lap();
        #endif
        

        //Add all the neighbors that were found
        auto neighbor_itr = params_.model_.get_neighbors(prev_kmer);
        for (auto next_kmer = neighbor_itr.first; 
             next_kmer != neighbor_itr.second; 
             next_kmer++) {

            prob = kmer_probs[*next_kmer];

            if (prob < evpr_thresh) {
                continue;
            }

            Base next_base = params_.model_.get_first_base(*next_kmer);

            #ifdef VERBOSE_TIME
            //neighbor_time += timer.lap();
            timer.reset();
            #endif

            Range next_range = fmi_.get_neighbor(prev_range, next_base);

            #ifdef VERBOSE_TIME
            fm_time += timer.lap();
            #endif

            child_found = 
                add_child(next_range,
                          prev_path, 
                          *next_kmer, 
                          prob, 
                          EventType::MATCH, 
                          child_found);

            //std::cout << "\t" << next_range.start_ << "-" 
            //          << next_range.end_ << " "
            //          << *next_kmer << "\n";

            #ifdef VERBOSE_TIME
            neighbor_time += timer.lap();
            #endif
        }

        #ifdef VERBOSE_TIME
        //neighbor_time += timer.lap();
        timer.reset();
        #endif


        if (child_found) {
            prev_path->update_consec_stays();
        } else {

            if (!prev_path->sa_checked_ && prev_path->should_report(prev_range, params_, false)) {
                Result r(cur_event_+1, params_.path_win_len_, prev_path->mean_prob());

                for (unsigned int s = prev_range.start_; s <= prev_range.end_; s++) {
                    r.set_ref_range(fmi_.sa(s), prev_path->match_len());
                    results.push_back(r);

                    seeds_out << label_ << "\t";
                    r.print(seeds_out);
                }
            }

            inactive_paths_.push_back(prev_path);
        }
        #ifdef VERBOSE_TIME
        un_ext_time += timer.lap();
        #endif
    }

    #ifdef VERBOSE_TIME
    time_out << std::setw(23) << stay_time 
             << std::setw(23) << fm_time
             << std::setw(23) << neighbor_time
             << std::setw(23) << un_ext_time
             << std::setw(23) << child_map_time_
             << std::setw(23) << child_add_time_
             << std::setw(23) << child_rpl_time_;
    #endif

    for (auto p = next_paths_.begin(); p != next_paths_.end(); p++) {

        PathBuffer *n = p->second;
        const Range &range = p->first;

        if (n->should_report(range, params_, true)) {
            Result r(cur_event_, params_.path_win_len_, n->mean_prob());
            n->sa_checked_ = true;

            for (unsigned int s = range.start_; s <= range.end_; s++) {
                r.set_ref_range(fmi_.sa(s), n->match_len());
                results.push_back(r);

                seeds_out << label_ << "\t";

                r.print(seeds_out);
            }
        }
    }

    //SA time
    #ifdef VERBOSE_TIME
    time_out << std::setw(23) << timer.lap();
    #endif

    //Find sources
    for (Kmer kmer = 0; kmer < params_.model_.kmer_count(); kmer++) {
        prob = kmer_probs[kmer];
        //std::cout << cur_event_ << " " << kmer << " " << params_.get_source_prob();
        if (prob >= params_.get_source_prob()) {
            Range next_range = kmer_ranges_[kmer];
            //std::cout << "\t" << next_range.start_ << "-" << next_range.end_ ;

            if (next_range.is_valid()) {
                //Will split sources if they intersect existing nodes
                add_sources(next_range, kmer, prob);             
            }
        }
        //std::cout << "\n";
    }
    
    //Source time
    #ifdef VERBOSE_TIME
    time_out << std::setw(23) << timer.lap();
    #endif

    prev_paths_.clear();
    prev_paths_.swap(next_paths_);

    return results;
}

size_t LeafAligner::add_sources(const Range &range, Kmer kmer, float prob) {

    //Find closest existing node
    auto lb = next_paths_.lower_bound(range);
    
    //Find range of existing alns intersecting the new aln
    auto start = lb;
    while (start != next_paths_.begin() 
           && std::prev(start)->first.intersects(range)){
        start--;
    }
    //if (!start->first.intersects(range)) {
    //    start++;
    //} 
    

    auto end = lb;
    while (end != next_paths_.end() && end->first.intersects(range)) {
        end++;
    }

    //No alns intersect new aln, just add it
    if (start == next_paths_.end() || (start == lb && end == lb)) {

        if (inactive_paths_.empty()) {
            next_paths_[range] = new PathBuffer(kmer, prob);
        } else {
            next_paths_[range] = inactive_paths_.back();
            next_paths_[range]->init_source(kmer, prob);
            inactive_paths_.pop_back();
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
        if (inactive_paths_.empty()) {
            next_paths_[*r] = new PathBuffer(kmer, prob);
        } else {
            next_paths_[*r] = inactive_paths_.back();
            next_paths_[*r]->init_source(kmer, prob);
            inactive_paths_.pop_back();
        }
    }

    return split_ranges.size();
}

bool LeafAligner::add_child(Range &range, 
                            PathBuffer *prev_path,
                            Kmer kmer,
                            float prob,
                            EventType type,
                            bool prev_is_sibling) {
    
    #ifdef VERBOSE_TIME
    Timer timer;
    timer.reset();
    #endif

    if (!range.is_valid()) {
        return prev_is_sibling;
    }

    
    //Find closest node >= the node being added
    auto lb = next_paths_.lower_bound(range);

    #ifdef VERBOSE_TIME
    child_map_time_ += timer.lap();
    #endif

    //std::cout << "a\n";

    //PathBuffer range hasn't been added yet
    if (lb == next_paths_.end() || !lb->first.same_range(range)) {

        PathBuffer *new_path;
        if (prev_is_sibling) {
            if (inactive_paths_.empty()) {
                //std::cout << "b\n";
                new_path = new PathBuffer(prev_path, kmer, prob, type);
            } else {
                //std::cout << "c\n";
                new_path = inactive_paths_.back();
                inactive_paths_.pop_back();
                new_path->init_from_sibling(prev_path, kmer, prob, type);
            }
        } else {
            //std::cout << "c/b\n";
            prev_path->make_child(kmer, prob, type);
            new_path = prev_path;
        }

        #ifdef VERBOSE_TIME
        child_add_time_ += timer.lap();
        #endif

        //Store it with it's range
        next_paths_.insert(lb, std::pair<Range, PathBuffer *>(range, new_path));

        #ifdef VERBOSE_TIME
        child_map_time_ += timer.lap();
        #endif

        //Created one node
        return true;
    }
    //std::cout << "old... \n";

    
    //PathBuffer associated with same range
    PathBuffer *dup_path = lb->second;

              
    //std::cout << "d " << prev_is_sibling << "\n";

    if (prev_is_sibling) {
        if (dup_path->better_than_sibling(prev_path, prob)) {
            dup_path->init_from_sibling(prev_path, kmer, prob, type);

            #ifdef VERBOSE_TIME
            child_rpl_time_ += timer.lap();
            #endif

            return true;
        }
    } else {
        if (dup_path->better_than_parent(prev_path, prob)) {
            inactive_paths_.push_back(lb->second);
            lb->second = prev_path;
            lb->second->make_child(kmer, prob, type);

            #ifdef VERBOSE_TIME
            child_rpl_time_ += timer.lap();
            #endif

            return true;
        }
    }

    #ifdef VERBOSE_TIME
    child_rpl_time_ += timer.lap();
    #endif

    return prev_is_sibling;
}

