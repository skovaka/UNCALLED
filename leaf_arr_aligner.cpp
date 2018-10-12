#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <stddef.h>
#include <utility>
#include "timer.hpp"
#include "fmi.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)
//#define DEBUG_PATHS
//#define DEBUG_RANGES

unsigned char LeafArrAligner::PathBuffer::MAX_WIN_LEN = 0, 
              LeafArrAligner::PathBuffer::TYPE_MASK = 0;

unsigned long LeafArrAligner::PathBuffer::TYPE_ADDS[EventType::NUM_TYPES];

LeafArrAligner::PathBuffer::PathBuffer()
    : length_(0),
      prob_sums_(new float[MAX_WIN_LEN+1]) {
}

LeafArrAligner::PathBuffer::PathBuffer(const PathBuffer &p) {
    std::memcpy(this, &p, sizeof(PathBuffer));
}

void LeafArrAligner::PathBuffer::free_buffers() {
    delete[] prob_sums_;
}

void LeafArrAligner::PathBuffer::make_source(Range &range, u16 kmer, float prob) {
    length_ = 1;
    consec_stays_ = 0;
    event_types_ = 0;
    win_prob_ = prob;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = false;

    //MATCH should always be index 0
    //all_type_counts_[EventType::MATCH] = 1;
    win_type_counts_[EventType::MATCH] = 1;
    for (unsigned char t = 1; t < EventType::NUM_TYPES; t++) {
        //all_type_counts_[t] = 0;
        win_type_counts_[t] = 0;
    }

    //TODO: don't write this here to speed up source loop
    prob_sums_[0] = 0;
    prob_sums_[1] = prob;
}


void LeafArrAligner::PathBuffer::make_child(PathBuffer &p, 
                                            Range &range,
                                            u16 kmer, 
                                            float prob, 
                                            EventType type) {
    length_ = p.length_ + 1;

    consec_stays_ = p.consec_stays_;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = p.sa_checked_;

    //std::memcpy(all_type_counts_, p.all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, p.win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));

    if (win_full()) {
        std::memcpy(prob_sums_, &(p.prob_sums_[1]), MAX_WIN_LEN * sizeof(float));
        prob_sums_[MAX_WIN_LEN] = prob_sums_[MAX_WIN_LEN-1] + prob;
        win_prob_ = (prob_sums_[MAX_WIN_LEN] - prob_sums_[0]) / MAX_WIN_LEN;
        win_type_counts_[p.type_tail()]--;
    } else {
        std::memcpy(prob_sums_, p.prob_sums_, (length_ + 1) * sizeof(float));
        prob_sums_[length_] = prob_sums_[length_-1] + prob;
        win_prob_ = (prob_sums_[length_] - prob_sums_[0]) / length_;
    }

    event_types_ = TYPE_ADDS[type] | (p.event_types_ >> TYPE_BITS);

    win_type_counts_[type]++;
    //all_type_counts_[type]++;

    if (type == EventType::STAY) {
        consec_stays_++;
    } else {
        consec_stays_ = 0;
    }
}

void LeafArrAligner::PathBuffer::invalidate() {
    length_ = 0;
}

bool LeafArrAligner::PathBuffer::is_valid() {
    return length_ > 0;
}

bool LeafArrAligner::PathBuffer::win_full() const {
    return length_ > MAX_WIN_LEN;
}

unsigned char LeafArrAligner::PathBuffer::win_len() const {
    return win_full() ? MAX_WIN_LEN : length_;
}

size_t LeafArrAligner::PathBuffer::event_len() {
    return length_;
}

size_t LeafArrAligner::PathBuffer::match_len() const {
    return win_type_counts_[EventType::MATCH];
}

float LeafArrAligner::PathBuffer::mean_prob() const {
    //return (prob_sums_[win_len()] - prob_sums_[0]) / win_len();
    return win_prob_;
}

bool LeafArrAligner::PathBuffer::better_than(const PathBuffer &p) {
    return p.mean_prob() > mean_prob();
}

unsigned char LeafArrAligner::PathBuffer::type_head() const {
    return (event_types_ >> (TYPE_BITS*(MAX_WIN_LEN-2))) & TYPE_MASK;
}

unsigned char LeafArrAligner::PathBuffer::type_tail() const {
    return event_types_ & TYPE_MASK;
}

bool LeafArrAligner::PathBuffer::should_report(const AlnParams &p,
                                                bool path_ended) {
    return (fm_range_.length() == 1 || 
                (path_ended &&
                 fm_range_.length() <= p.max_rep_copy_ &&
                 match_len() >= p.min_rep_len_)) &&

           length_ >= p.path_win_len_ &&
           (path_ended || type_head() == EventType::MATCH) &&
           (path_ended || win_type_counts_[EventType::STAY] <= p.max_stay_frac_ * p.path_win_len_) &&
          win_prob_ >= p.window_prob_;
}

void LeafArrAligner::PathBuffer::replace(const PathBuffer &p) {
    length_ = p.length_;
    win_prob_ = p.win_prob_;
    consec_stays_ = p.consec_stays_;
    fm_range_ = p.fm_range_;
    kmer_ = p.kmer_;
    sa_checked_ = p.sa_checked_;

    std::memcpy(prob_sums_, p.prob_sums_, (MAX_WIN_LEN+1) * sizeof(float));
    //std::memcpy(all_type_counts_, p.all_type_counts_, EventType::NUM_TYPES * sizeof(unsigned short));
    std::memcpy(win_type_counts_, p.win_type_counts_, EventType::NUM_TYPES * sizeof(unsigned char));
}

bool operator< (const LeafArrAligner::PathBuffer &p1, 
                const LeafArrAligner::PathBuffer &p2) {
    return p1.fm_range_ < p2.fm_range_ ||
           (p1.fm_range_ == p2.fm_range_ && 
            p1.win_prob_ < p2.win_prob_);
}

LeafArrAligner::LeafArrAligner(const FMI &fmi, 
                     const AlnParams &ap)
    : fmi_(fmi),
      params_(ap) {

    PathBuffer::MAX_WIN_LEN = params_.path_win_len_;

    for (unsigned long t = 0; t < EventType::NUM_TYPES; t++) {
        PathBuffer::TYPE_ADDS[t] = t << ((PathBuffer::MAX_WIN_LEN-2)*TYPE_BITS);
    }
    PathBuffer::TYPE_MASK = (unsigned char) ((1 << TYPE_BITS) - 1);

    kmer_ranges_ = new Range[params_.model_.kmer_count()];
    for (u16 k = 0; k < params_.model_.kmer_count(); k++) {
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


    #ifdef VERBOSE_TIME
    loop1_time_ = fmrs_time_ = fmsa_time_ = sort_time_ = loop2_time_ = fullsource_time_ = 0;
    #endif
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
    cur_event_ = 0;
    #ifdef VERBOSE_TIME
    loop1_time_ = fmrs_time_ = fmsa_time_ = sort_time_ = loop2_time_ = fullsource_time_ = 0;
    #endif
}

void LeafArrAligner::reset() {
    prev_size_ = 0;
}

std::vector<Result> LeafArrAligner::add_event(std::vector<float> kmer_probs, 
                                           std::ostream &seeds_out,
                                           std::ostream &time_out) {


    float prob;

    Range prev_range;
    //PathBuffer *prev_path;
    u16 prev_kmer;
    
    bool child_found;
    float evpr_thresh;

    std::vector<Result> results;

    #ifdef VERBOSE_TIME
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

        //Get probability for stay neighbor
        prob = kmer_probs[prev_kmer];

        if (prev_path.consec_stays_ < params_.max_consec_stay_ && 
            prob >= evpr_thresh) {

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
        for (u8 i = 0; i < 4; i++) {
            u16 next_kmer = params_.model_.get_neighbor(prev_kmer, i);

            prob = kmer_probs[next_kmer];

            if (prob < evpr_thresh) {
                continue;
            }

            //u8 next_base = params_.model_.get_first_base(*next_kmer);
            u8 next_base = params_.model_.get_last_base(next_kmer);

            #ifdef VERBOSE_TIME
            loop1_time_ += timer.lap();
            #endif

            Range next_range = fmi_.get_neighbor(prev_range, next_base);

            #ifdef VERBOSE_TIME
            fmrs_time_ += timer.lap();
            #endif

            if (!next_range.is_valid()) {
                continue;
            }

            next_path->make_child(prev_path, 
                                  next_range,
                                  next_kmer, 
                                  prob, 
                                  EventType::MATCH);

            child_found = true;

            if (++next_path == next_paths_.end()) {
                break;
            }
        }

        if (!child_found && !prev_path.sa_checked_) {
            #ifdef VERBOSE_TIME
            loop1_time_ += timer.lap();
            #endif

            check_alignments(prev_path, results, true, seeds_out);

            #ifdef VERBOSE_TIME
            fmsa_time_ += timer.lap();
            #endif
        }

        if (next_path == next_paths_.end()) {
            break;
        }
    }

    #ifdef VERBOSE_TIME
    loop1_time_ += timer.lap();
    #endif

    if (next_path != next_paths_.begin()) {

        size_t next_size = next_path - next_paths_.begin();

        std::sort(next_paths_.begin(), next_path);

        #ifdef VERBOSE_TIME
        sort_time_ += timer.lap();
        #endif

        u16 source_kmer;
        prev_kmer = params_.model_.kmer_count(); //TODO: won't work for k=4

        Range unchecked_range, source_range;

        for (u32 i = 0; i < next_size; i++) {
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

            //Range next_range = next_paths_[i].fm_range_;

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
                
                source_range = unchecked_range;
                
                //Between this and next path ranges
                if (i < next_size - 1 && source_kmer == next_paths_[i+1].kmer_) {

                    source_range.end_ = next_paths_[i+1].fm_range_.start_ - 1;

                    if (unchecked_range.start_ <= next_paths_[i+1].fm_range_.end_) {
                        unchecked_range.start_ = next_paths_[i+1].fm_range_.end_ + 1;
                    }
                }

                //Add it if it's a real range
                if (source_range.is_valid()) {
                    next_path->make_source(source_range,
                                           source_kmer,
                                           prob);
                    next_path++;
                }
            }

            #ifdef VERBOSE_TIME
            loop2_time_ += timer.lap();
            #endif

            check_alignments(next_paths_[i], results, false, seeds_out);

            #ifdef VERBOSE_TIME
            fmsa_time_ += timer.lap();
            #endif
        }
    }

    #ifdef VERBOSE_TIME
    loop2_time_ += timer.lap();
    #endif
    
    for (u16 kmer = 0; 
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

    #ifdef VERBOSE_TIME
    fullsource_time_ += timer.lap();
    #endif


    //std::cout << ((prev_size_ + (next_path - next_paths_.begin())) / 2) << "\t" << timer.get() << "\n";
    prev_size_ = next_path - next_paths_.begin();
    prev_paths_.swap(next_paths_);

    //Update event index
    cur_event_++;

    return results;
}

void LeafArrAligner::check_alignments(PathBuffer &p, 
                                      std::vector<Result> &results, 
                                      bool path_ended,
                                      std::ostream &seeds_out) {

    if (p.should_report(params_, path_ended)) {
        //std::cout << " pass\n";

        Result r(cur_event_ - path_ended, params_.path_win_len_, p.win_prob_);

        p.sa_checked_ = true;

        for (u64 s = p.fm_range_.start_; s <= p.fm_range_.end_; s++) {

            //Reverse the reference coords so they both go L->R
            u64 rev_en = fmi_.size() - fmi_.sa(s) + 1;
            r.set_ref_range(rev_en, p.match_len());

            //r.set_ref_range(fmi_.sa(s), p.match_len());
            
            results.push_back(r);

            //seeds_out << label_ << "\t";
            //r.print(seeds_out);
        }
    } else {
        //std::cout << " fail\n";
    }
}

