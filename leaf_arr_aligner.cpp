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
#include "leaf_arr_aligner.hpp"
#include "fmi.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)
//#define DEBUG_PATHS
//#define DEBUG_RANGES

#define MAX_CHILDREN 5

unsigned char LeafArrAligner::PathBuffer::PROB_WIN_LEN = 0, 
              LeafArrAligner::PathBuffer::TYPE_WIN_LEN = 0;

LeafArrAligner::PathBuffer::PathBuffer()
    : length_(0),
      prob_sums_(new float[PROB_WIN_LEN]),
      event_types_(new EventType[TYPE_WIN_LEN]) {
}

LeafArrAligner::PathBuffer::PathBuffer(const PathBuffer &p) :
    length_ (p.length_),
    consec_stays_ (p.consec_stays_),
    prhd_   (p.prhd_),
    prtl_   (p.prtl_),
    prlen_  (p.prlen_),
    tyhd_   (p.tyhd_),
    tytl_   (p.tytl_),
    tylen_  (p.tylen_),
    win_prob_ (p.win_prob_),
    prob_sums_ (p.prob_sums_),
    event_types_ (p.event_types_),
    fm_range_ (p.fm_range_),
    kmer_ (p.kmer_),
    sa_checked_ (p.sa_checked_) {

    //all_type_counts_ = p.all_type_counts_;
    //win_type_counts_ = p.win_type_counts_;

    //std::memcpy(prob_sums_, p.prob_sums_, PROB_WIN_LEN * sizeof(float));
    //std::memcpy(event_types_, p.event_types_, TYPE_WIN_LEN * sizeof(EventType));
    std::memcpy(all_type_counts_, p.all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, p.win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));
}

void LeafArrAligner::PathBuffer::free_buffers() {
    delete[] prob_sums_;
    delete[] event_types_;
}
//LeafArrAligner::PathBuffer::~PathBuffer() {
//    std::cout << "DESTRUCT\n";
//}

void LeafArrAligner::PathBuffer::make_source(Range &range, Kmer kmer, float prob) {
    length_ = 1;
    consec_stays_ = 0;
    prtl_ = 0;
    prhd_ = 1;
    prlen_ = 1;
    tytl_ = 0;
    tyhd_ = 0;
    tylen_ = 1;
    win_prob_ = prob;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = false;

    //MATCH should always be index 0
    all_type_counts_[EventType::MATCH] = 1;
    win_type_counts_[EventType::MATCH] = 1;
    for (unsigned char t = 1; t < EventType::NUM_TYPES; t++) {
        all_type_counts_[t] = 0;
        win_type_counts_[t] = 0;
    }

    //TODO: don't write this here to speed up source loop
    prob_sums_[prtl_] = 0;
    prob_sums_[prhd_] = prob;
    event_types_[tyhd_] = EventType::MATCH;
}

/*
void LeafArrAligner::PathBuffer::init_from_sibling(PathBuffer *a, Kmer kmer, 
                                               float prob, EventType type) {
    length_ = a->length_;
    prtl_ = a->prtl_;
    prhd_ = a->prhd_;
    prlen_ = a->prlen_;
    tytl_ = a->tytl_;
    tyhd_ = a->tyhd_;
    tylen_ = a->tylen_;
    kmer_ = kmer;
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
*/

void LeafArrAligner::PathBuffer::make_child(PathBuffer &p, 
                                            Range &range,
                                            Kmer kmer, 
                                            float prob, 
                                            EventType type) {
    length_ = p.length_ + 1;

    //TODO: float initing some of these below
    prtl_ = p.prtl_;
    prhd_ = p.prhd_;
    prlen_ = p.prlen_;
    tytl_ = p.tytl_;
    tyhd_ = p.tyhd_;
    tylen_ = p.tylen_;

    consec_stays_ = p.consec_stays_;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = p.sa_checked_;

    //std::cout << prob_sums_ << " " << p.prob_sums_ << " " << length_ << "\n";
    std::memcpy(prob_sums_, p.prob_sums_, PROB_WIN_LEN * sizeof(float));
    std::memcpy(event_types_, p.event_types_, TYPE_WIN_LEN * sizeof(EventType));
    std::memcpy(all_type_counts_, p.all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, p.win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));

    if (tylen_ < TYPE_WIN_LEN) { //TODO: maybe off by 1
        tyhd_++;
        tylen_++;
    } else {
        tyhd_ = tytl_;
        tytl_ = (tytl_ + 1) % TYPE_WIN_LEN;
        win_type_counts_[event_types_[tyhd_]]--;

        win_type_counts_[event_types_[tytl_]]--;
        event_types_[tytl_] = EventType::MATCH;
        win_type_counts_[EventType::MATCH]++;
    }

    event_types_[tyhd_] = type;
    win_type_counts_[type]++;
    all_type_counts_[type]++;

    if (type == EventType::STAY) {
        consec_stays_++;
    } else {
        consec_stays_ = 0;
    }

    if (prlen_ < PROB_WIN_LEN - 1) { //TODO: maybe off by 1
        prob_sums_[prhd_+1] = prob_sums_[prhd_] + prob;
        prhd_++;
        prlen_++;
    } else {
        prob_sums_[prtl_] = prob_sums_[prhd_] + prob;
        prhd_ = prtl_;
        prtl_ = (prtl_ + 1) % PROB_WIN_LEN;
    }

    win_prob_ = mean_prob();
}

void LeafArrAligner::PathBuffer::invalidate() {
    length_ = 0;
}

bool LeafArrAligner::PathBuffer::is_valid() {
    return length_ > 0;
}

size_t LeafArrAligner::PathBuffer::event_len() {
    return length_;
}

size_t LeafArrAligner::PathBuffer::match_len() {
    return win_type_counts_[EventType::MATCH];
}

float LeafArrAligner::PathBuffer::mean_prob() const {
    return (prob_sums_[prhd_] - prob_sums_[prtl_]) / prlen_;
}

bool LeafArrAligner::PathBuffer::better_than(const PathBuffer &p) {
    return p.mean_prob() > mean_prob();
}

/*
bool LeafArrAligner::PathBuffer::better_than_sibling(const PathBuffer *a, float prob) {
    float x = a->prob_sums_[a->prhd_ > 0 ? a->prhd_ - 1 : PROB_WIN_LEN-1],
           y = a->prob_sums_[a->prtl_];

    float replace_prob = (x - y + prob) / a->prlen_;

    return replace_prob > mean_prob();
}

bool LeafArrAligner::PathBuffer::better_than_parent(const PathBuffer *a, float prob) {
    float x = a->prob_sums_[a->prhd_], 
           y = a->prob_sums_[(a->prlen_ == PROB_WIN_LEN-1 && a->prtl_ < PROB_WIN_LEN-1) ? a->prtl_ + 1 : 0];
    unsigned int len = a->prlen_ + (a->prlen_ < PROB_WIN_LEN-1);
    float replace_prob = (x - y + prob) / len;
    //std::cout << (int) a->prtl_ << "\t" << (int) a->prhd_ << "\n";

    return replace_prob > mean_prob();
}
*/

float LeafArrAligner::PathBuffer::next_mean_prob() {
    float x = prob_sums_[prhd_], 
           y = prob_sums_[(prlen_ == PROB_WIN_LEN-1 && prtl_ < PROB_WIN_LEN-1) ? prtl_ + 1 : 0];
    unsigned int len = prlen_ - (prlen_ == PROB_WIN_LEN-1);
    return (x - y) / len;
}

bool LeafArrAligner::PathBuffer::should_report(const AlnParams &p,
                                                bool path_ended) {
    //if (length_ >= 22) {
    //    std::cout << length_ << " " 
    //              << fm_range_.length() << " "
    //              << event_types_[tyhd_] << " "
    //              << mean_prob() << " "
    //              << (int) win_type_counts_[EventType::STAY] << "\n";
    //    std::cout.flush();
    //}
    //
    return (fm_range_.length() == 1 || 
                (path_ended &&
                 fm_range_.length() <= p.max_rep_copy_ &&
                 match_len() >= p.min_rep_len_)) &&

           length_ >= p.path_win_len_ &&
           event_types_[tyhd_] == EventType::MATCH &&
           (path_ended || win_type_counts_[EventType::STAY] <= p.max_stay_frac_ * tylen_) &&
           win_prob_ >= p.window_prob_;
}

void LeafArrAligner::PathBuffer::replace(const PathBuffer &p) {
    length_ = p.length_;
    prtl_ = p.prtl_;
    prhd_ = p.prhd_;
    prlen_ = p.prlen_;
    tytl_ = p.tytl_;
    tyhd_ = p.tyhd_;
    tylen_ = p.tylen_;
    win_prob_ = p.win_prob_;
    consec_stays_ = p.consec_stays_;
    fm_range_ = p.fm_range_;
    kmer_ = p.kmer_;
    sa_checked_ = p.sa_checked_;

    std::memcpy(prob_sums_, p.prob_sums_, PROB_WIN_LEN * sizeof(float));
    std::memcpy(event_types_, p.event_types_, TYPE_WIN_LEN * sizeof(EventType));
    std::memcpy(all_type_counts_, p.all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, p.win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));
}

bool operator< (const LeafArrAligner::PathBuffer &p1, 
                const LeafArrAligner::PathBuffer &p2) {
    return p1.fm_range_ < p2.fm_range_ ||
           (p1.fm_range_ == p2.fm_range_ && 
            p1.win_prob_ < p2.win_prob_);
}

LeafArrAligner::LeafArrAligner(const FMI &fmi, 
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
        //std::cout << k << " " << r.start_ << "-" << r.end_ << "\n";
        kmer_ranges_[k] = r;
    }

    prev_paths_ = std::vector<PathBuffer>(params_.max_paths_);
    next_paths_ = std::vector<PathBuffer>(params_.max_paths_);
    sources_added_ = std::vector<bool>(params_.model_.kmer_count(), false);
    prev_size_ = 0;
}

LeafArrAligner::~LeafArrAligner() {
    delete[] kmer_ranges_;
    for (size_t i = 0; i < next_paths_.size(); i++) {
        next_paths_[i].free_buffers();
        prev_paths_[i].free_buffers();
    }
    reset();
}

void LeafArrAligner::new_read(size_t read_len) {
    reset();
    cur_event_ = read_len;
}

void LeafArrAligner::reset() {
    prev_size_ = 0;
}

std::vector<Result> LeafArrAligner::add_event(float *kmer_probs, 
                                           std::ostream &seeds_out,
                                           std::ostream &time_out) {

    //Update event index
    cur_event_--;

    float prob;

    Range prev_range;
    //PathBuffer *prev_path;
    Kmer prev_kmer;
    
    bool child_found;
    float evpr_thresh;

    std::vector<Result> results;

    #ifdef VERBOSE_TIME
    //float stay_time = 0, 
    //       fm_time = 0, 
    //       neighbor_time = 0, 
    //       un_ext_time = 0;
    float first_loop_time = 0,
           sort_time = 0,
           second_loop_time = 0;
    Timer timer;
    #endif

    auto next_path = next_paths_.begin();

    //Find neighbors of previous nodes
    for (size_t pi = 0; pi < prev_size_; pi++) {
        if (!prev_paths_[pi].is_valid()) {
            continue;
        }


        child_found = false;

        PathBuffer &prev_path = prev_paths_[pi];
        Range &prev_range = prev_path.fm_range_;
        prev_kmer = prev_path.kmer_;

        evpr_thresh = params_.get_prob_thresh(prev_range.length());

        //std::cout << cur_event_ << "\t" 
        //          << prev_path.length_ << "\t"
        //          << prev_path.mean_prob() << "\t"
        //          << prev_kmer << "\t"
        //          << prev_path.event_types_[prev_path.tyhd_] << "\t"
        //          << prev_path.fm_range_.start_ << "-"
        //          << prev_path.fm_range_.end_ << "\n";


        //#ifdef VERBOSE_TIME
        //timer.reset();
        //#endif

        //Get probability for stay neighbor
        prob = kmer_probs[prev_kmer];

        if (prev_path.consec_stays_ < params_.max_consec_stay_ && 
            prob >= evpr_thresh) {

            //std::cout << (next_paths_.end() - next_path) << "\n";
            next_path->make_child(prev_path, 
                           prev_range,
                           prev_kmer, 
                           prob, 
                           EventType::STAY);

            child_found = true;

            if (++next_path == next_paths_.end()) {
                break;
            }
        }

        //#ifdef VERBOSE_TIME
        //stay_time += timer.lap();
        //#endif
        
        //Add all the neighbors
        auto neighbor_itr = params_.model_.get_neighbors(prev_kmer);
        for (auto next_kmer = neighbor_itr.first; 
             next_kmer != neighbor_itr.second; 
             next_kmer++) {

            prob = kmer_probs[*next_kmer];

            if (prob < evpr_thresh) {
                continue;
            }

            Base next_base = params_.model_.get_first_base(*next_kmer);

            //#ifdef VERBOSE_TIME
            ////neighbor_time += timer.lap();
            //timer.reset();
            //#endif

            Range next_range = fmi_.get_neighbor(prev_range, next_base);

            //#ifdef VERBOSE_TIME
            //fm_time += timer.lap();
            //#endif

            if (!next_range.is_valid()) {
                continue;
            }

            //std::cout << (next_path == next_paths_.end()) << "\n";

            next_path->make_child(prev_path, 
                                  next_range,
                                  *next_kmer, 
                                  prob, 
                                  EventType::MATCH);

            child_found = true;

            if (++next_path == next_paths_.end()) {
                break;
            }

            //#ifdef VERBOSE_TIME
            //neighbor_time += timer.lap();
            //#endif
        }


        //#ifdef VERBOSE_TIME
        //neighbor_time += timer.lap();
        //#endif


        if (!child_found && !prev_path.sa_checked_) {
            check_alignments(prev_path, results, true);
        }

        //#ifdef VERBOSE_TIME
        //un_ext_time += timer.lap();
        //#endif

        if (next_path == next_paths_.end()) {
            break;
        }
    }

    #ifdef VERBOSE_TIME
    time_out << std::setw(23) << timer.lap() << "\t";
    #endif

    if (next_path != next_paths_.begin()) {
        size_t next_size = next_path - next_paths_.begin();

        std::sort(next_paths_.begin(), next_path);

        #ifdef VERBOSE_TIME
        time_out << std::setw(23) << timer.lap() << "\t";
        #endif

        //std::vector<PathBuffer *> paths_ptrs(next_size);
        //for (size_t i = 0; i < next_size; i++) {
        //    paths_ptrs[i] = &(next_paths_[i]);
        //}

        //Timer t;
        //std::sort(paths_ptrs.begin(), paths_ptrs.end(),
        //          [](const PathBuffer *p1, const PathBuffer *p2) -> bool {
        //              return *p1 < *p2;
        //          });
        //std::cout << next_size << "\t" << t.lap() << "\t";
        //std::cout << t.lap() << "\t" << sizeof(next_paths_[0]) << "\n";


        Kmer source_kmer;
        prev_kmer = params_.model_.kmer_count(); //TODO: won't work for k=4

        Range unchecked_range, source_range;

        for (size_t i = 0; i < next_size; i++) {
            source_kmer = next_paths_[i].kmer_;
            prob = kmer_probs[source_kmer];

            //Add source for beginning of kmer range
            if (source_kmer != prev_kmer &&
                next_path != next_paths_.end() &&
                prob >= params_.get_source_prob()) {

                sources_added_[source_kmer] = true;

                //unchecked_range = kmer_ranges_[source_kmer];

                source_range = Range(kmer_ranges_[source_kmer].start_,
                                     next_paths_[i].fm_range_.start_ - 1);

                if (source_range.is_valid()) {
                    next_path->make_source(source_range,
                                           source_kmer,
                                           prob);
                    next_path++;
                }                                    

                unchecked_range = Range(next_paths_[i].fm_range_.end_ + 1,
                                        kmer_ranges_[source_kmer].end_);
            }

            prev_kmer = source_kmer;

            //Remove paths with duplicate ranges
            //Best path will be listed last
            if (i < next_size - 1 && next_paths_[i].fm_range_ == next_paths_[i+1].fm_range_) {
                next_paths_[i].invalidate();
                continue;
            }

            //Start source after current path
            //TODO: check if theres space for a source here, instead of after extra work?
            if (next_path != next_paths_.end() &&
                prob >= params_.get_source_prob()) {
                
                //if (unchecked_range.start_ <= next_paths_[i].fm_range_.end_) {
                //    unchecked_range.start_ = next_paths_[i].fm_range_.end_ + 1; 
                //}

                source_range = unchecked_range;
                
                //Between this and next path ranges
                if (i < next_size - 1 && source_kmer == next_paths_[i+1].kmer_) {

                    source_range.end_ = next_paths_[i+1].fm_range_.start_ - 1;

                    if (unchecked_range.start_ <= next_paths_[i+1].fm_range_.end_) {
                        unchecked_range.start_ = next_paths_[i+1].fm_range_.end_ + 1;
                    }

                //Between this path and end of kmer
                }// else {
                //    source_range.end_ = kmer_ranges_[source_kmer].end_;
                //}

                //Add it if it's a real range
                if (source_range.is_valid()) {
                    next_path->make_source(source_range,
                                           source_kmer,
                                           prob);
                    next_path++;
                }
            }

            check_alignments(next_paths_[i], results, false);
        }
    }
    
    #ifdef VERBOSE_TIME
    else {
    time_out << std::setw(23) << timer.lap() << "\t";
    }
    #endif

    #ifdef VERBOSE_TIME
    time_out << std::setw(23) << timer.lap() << "\t";
    #endif

    for (Kmer kmer = 0; 
         kmer < params_.model_.kmer_count() && 
            next_path != next_paths_.end(); 
         kmer++) {

        prob = kmer_probs[kmer];
        Range next_range = kmer_ranges_[kmer];

        if (!sources_added_[kmer] && 
            prob >= params_.get_source_prob() &&
            next_path != next_paths_.end() &&
            next_range.is_valid()) {

            //TODO: don't write to prob buffer here to speed up source loop
            next_path->make_source(next_range, kmer, prob);             
            next_path++;

        } else {
            sources_added_[kmer] = false;
        }
    }

    //SA time
    #ifdef VERBOSE_TIME
    time_out << std::setw(23) << timer.lap();
    #endif


    prev_size_ = next_path - next_paths_.begin();
    //if (prev_size_ == next_paths_.size()) {
    //    std::cout << "FULL\n";
    //}
    //std::cout << prev_size_ << "\n";
    prev_paths_.swap(next_paths_);

    return results;
}

void LeafArrAligner::check_alignments(PathBuffer &p, 
                                      std::vector<Result> &results, 
                                      bool path_ended) {

    if (p.should_report(params_, path_ended)) {
        Result r(cur_event_ + path_ended, params_.path_win_len_, p.win_prob_);
        p.sa_checked_ = true;

        for (unsigned int s = p.fm_range_.start_; s <= p.fm_range_.end_; s++) {
            r.set_ref_range(fmi_.sa(s), p.match_len());
            results.push_back(r);

            //seeds_out << label_ << "\t";

            //r.print(seeds_out);
        }
    }
}

