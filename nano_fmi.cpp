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
#include "timer.h"
#include "nano_fmi.hpp"
#include "boost/math/distributions/students_t.hpp"

//Reads a model directly from a file and creates the FM index from the given reference
NanoFMI::NanoFMI(KmerModel &model, std::vector<mer_id> &mer_seq, int tally_dist) {

    model_ = &model;
    mer_seq_ = &mer_seq;
    tally_dist_ = tally_dist;

    //For outputting time
    Timer timer;

    //Init suffix array
    //Not using suffix_ar instance var speeds up sorting significantly
    std::vector<unsigned int> suffix_ar(mer_seq.size());
    for (unsigned int i = 0; i < suffix_ar.size(); i++) 
        suffix_ar[i] = i;

    std::cerr << "SA init time: " << timer.lap() << "\n";

    //Create the suffix array
    std::sort(suffix_ar.begin(), suffix_ar.end(), *this);
    suffix_ar_.swap(suffix_ar);

    std::cerr << "SA sort time: " << timer.lap() << "\n";

    //Allocate space for other datastructures
    bwt_.resize(mer_seq.size());
    mer_f_starts_.resize(model.kmer_count());
    mer_counts_.resize(model.kmer_count());
    mer_tally_.resize(model.kmer_count());

    for (mer_id i = 0; i < model.kmer_count(); i++)
        mer_tally_[i].resize((mer_seq.size() / tally_dist_) + 1, -1);
    
    std::cerr << "FM init time: " << timer.lap() << "\n";

    int tally_mod = tally_dist_;
    
    //Single pass to generate BWT and other datastructures
    for (unsigned int i = 0; i < suffix_ar_.size(); i++) {
        
        //Fill in BWT
        if (suffix_ar_[i] > 0)
            bwt_[i] = mer_seq[suffix_ar_[i]-1];
        else
            bwt_[i] = mer_seq[suffix_ar_[suffix_ar_.size()-1]];

        //Update 6-mer counts
        mer_counts_[bwt_[i]]++;
        
        //Update tally array
        if (tally_mod == tally_dist_) {
            for (mer_id j = 0; j < model.kmer_count(); j++)
                mer_tally_[j][i / tally_dist_] = mer_counts_[j];
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    std::cerr << "FM build time: " << timer.lap() << "\n";
    
    //Compute start locations for F array
    for (mer_id i = 0; i < model.kmer_count(); i++) {
        mer_f_starts_[i] = 1;
        for (mer_id j = 0; j < model.kmer_count(); j++)
            if (model_->compare_kmers(i, j) > 0)
                mer_f_starts_[i] += mer_counts_[j];
    }
    
    //Fill in last entry in tally array if needed
    if (mer_seq.size() % tally_dist_ == 0)
        for (mer_id i = 0; i < model.kmer_count(); i++)
            mer_tally_[i][mer_tally_[i].size()-1] = mer_counts_[i];

}

//Returns true if the suffix of *mer_seq_tmp starting at rot1 is less than that
//starting at rot2. Used to build suffix array.
bool NanoFMI::operator() (unsigned int rot1, unsigned int rot2) {
    
    int c;
    for (unsigned int i = 0; i < mer_seq_->size(); i++) {
        c = model_->compare_kmers(mer_seq_->at(rot1+i), mer_seq_->at(rot2+i));
        
        if (c == 0)
            continue;

        if (c < 0)
            return true;

       return false;
    }

    return false;
}


float NanoFMI::get_stay_prob(Event e1, Event e2) {
    double var1 = e1.stdv*e1.stdv, var2 = e2.stdv*e2.stdv;

    double t = (e1.mean - e2.mean) / sqrt(var1/e1.length + var2/e2.length);

    int df = pow(var1/e1.length + var2/e2.length, 2) 
               / (pow(var1/e1.length, 2) / (e1.length-1) 
                    + pow(var2/e2.length, 2) / (e2.length-1));

    boost::math::students_t dist(df);
    double q = boost::math::cdf(boost::math::complement(dist, fabs(t)));

    return q;
}

//Returns the distance from a BWT index to the nearest tally array checkpoint
int NanoFMI::tally_cp_dist(int i) {
    int cp = (i / tally_dist_)*tally_dist_; //Closest checkpoint < i

    //Check if checkpoint after i is closer
    if (i - cp > (cp + tally_dist_) - i && cp + (unsigned) tally_dist_ < bwt_.size())
        cp += tally_dist_;

    return cp - i;
}

//Returns the number of occurences of the given k-mer in the BWT up to and
//including the given index
int NanoFMI::get_tally(mer_id c, int i) {
    if (i < 0)
        return -1;
    int cp_dist = tally_cp_dist(i);
    int tally = mer_tally_[c][(i + cp_dist) / tally_dist_];

    if (cp_dist > 0) {
        for (int j = i+1; j <= i + cp_dist; j++)
            if (bwt_[j] == c)
                tally--;

    } else if (cp_dist < 0) {
        for (int j = i; j > i + cp_dist; j--)
            if (bwt_[j] == c)
                tally++;
    }

    return tally;
}

//Aligns the vector of events to the reference using LF mapping
//Returns the number of exact alignments
std::vector<NanoFMI::Result> NanoFMI::lf_map(
                    std::vector<Event> &events, int map_start, 
                    int seed_len, NormParams norm_params) {

    //float EVENT_THRESH = 0.00025,
    //      SEED_THRESH  = 0.05,
    //      STAY_THRESH  = 0.01;

    float EVENT_THRESH = -9.2103,
          SEED_THRESH = -3.75, //-5.298,
          STAY_THRESH = -5.298;


    //Match the first event

    double **event_probs = new double*[map_start + 1];//[model_->kmer_count()];

    std::cerr << "Getting probs\n";

    for (int e = map_start; e >= 0; e--) {
        event_probs[e] = new double[model_->kmer_count()];
        for (mer_id i = 0; i < model_->kmer_count(); i++) {
            event_probs[e][i] = model_->event_match_prob(events[e], i, norm_params);
        }
    }
    std::cerr << "Got probs\n";

    std::list< std::set<Range> > traversed_ranges;
    //traversed_ranges.push_front(std::set<Range>());
    //std::set<Range> to_traverse;

    std::vector<Result> results;

    //double prob;

    int i = map_start;
    while (i >= 0) {

        //std::cout << "loop " << i << "\n";

        std::set<Range> next_ranges, seed_starts;
        
        bool stay = false;

        //stay = get_stay_prob(events[i], events[i+1]) >= STAY_THRESH;
        stay = true;

        //auto prev_ranges = traversed_ranges.begin(); //Empty on first iteration
        if (!traversed_ranges.empty()) {
            for (auto pr = traversed_ranges.front().begin(); 
                 pr != traversed_ranges.front().end(); pr++) {

                if (stay && event_probs[i][pr->k_id_] >= EVENT_THRESH) {
                    Range nr(*pr, event_probs[i][pr->k_id_]);

                    update_ranges(next_ranges, nr, false);
                }
                
                auto neighbors = model_->get_neighbors(pr->k_id_);

                for (auto n = neighbors.first; n != neighbors.second; n++) {

                    if(event_probs[i][*n] >= EVENT_THRESH) {

                        Range nr(*pr, *n, event_probs[i][*n]);
                        //std::cout << nr.match_len_ << "\n";
                        if (update_ranges(next_ranges, nr, false) > 0 
                                 && nr.match_len_ == seed_len ) {
                            update_ranges(seed_starts, nr, false);
                        }
                    }
                }
            }
        }

        //std::cout << next_ranges.size() << ", ";

        //FIRST MATCH
        for (mer_id k_id = 0; k_id < model_->kmer_count(); k_id++) {
            if (event_probs[i][k_id] >= EVENT_THRESH) {
                Range nr(*this, k_id, i, event_probs[i][k_id]);
                update_ranges(next_ranges, nr, true);
                //std::cout << update_ranges(next_ranges, nr, true) << "\n";
            }
        }

       // std::cout << next_ranges.size() << "\n";

        for (auto aln_st = seed_starts.begin(); aln_st != seed_starts.end(); aln_st++) {

            auto end_ranges = traversed_ranges.begin();
            int total_length = 1;
            bool back_stays = false;
            const Range *aln_en = &(*aln_st),
                        *aln_cut_end = NULL; //Marks

            while (aln_en->parent_ != NULL) {
                if (!back_stays) //Last 
                    aln_cut_end = aln_en;
                back_stays = aln_en->stay_;
                aln_en = aln_en->parent_;
                total_length++;
                end_ranges++;
            }
                
            if (aln_st->prob_ / total_length >= SEED_THRESH) {
                for (int s = aln_st->start_; s <= aln_st->end_; s++) {
                    Result r;
                    r.qry_start = aln_st->event_;// seed_end - total_length + 1;
                    r.qry_end = aln_en->event_;
                    r.ref_start = suffix_ar_[s];
                    r.ref_end = suffix_ar_[s] + aln_st->match_len_ - 1;
                    r.prob = aln_st->prob_ / aln_st->total_len_;
                    results.push_back(r);

                    std::cout << "rev\t" 
                              << r.qry_start << "-" << r.qry_end << "\t"
                              << r.ref_start << "-" << r.ref_end << "\t"
                              << r.prob << "\n";//" " << total_length << " "
                              //<< aln_st->match_len_ << " " <<  aln_st->event_ << "\n";
                }
            }

            //double last_prob = event_probs[aln_en->event_][aln_en->k_id_];
            //double last_prob = aln_en->prob_;
            double last_prob = aln_cut_end->prob_;

            aln_st->match_len_--;
            aln_st->prob_ -= last_prob;

            int next_len = aln_st->match_len_ - 1;

            const Range *aln_evt = &(*aln_st);
            while ((aln_evt = aln_evt->parent_) != NULL) {

                if (aln_evt->match_len_ != next_len) {
                    aln_evt->prob_ -= last_prob;
                    aln_evt->match_len_ = next_len;
                    aln_evt->total_len_--;
                }

                if (aln_evt->parent_ == aln_cut_end) {
                    aln_evt->parent_ = NULL;
                    aln_cut_end->child_count_--;
                }

                else if (aln_evt->stay_)
                    next_len = aln_evt->match_len_;
                else
                    next_len = aln_evt->match_len_-1;
            }

            
            aln_evt = aln_cut_end;
            while (aln_evt != NULL && aln_evt->child_count_ == 0) {
                if ((aln_evt = aln_evt->parent_) != NULL) 
                    aln_evt->child_count_--;
            }
        }

        //std::cout << traversed_ranges.size() << " " << next_ranges.size() << "\n";

        traversed_ranges.push_front(next_ranges);

        i--;
    }

    delete event_probs;

    return results;
}

int NanoFMI::update_ranges(std::set<NanoFMI::Range> &next_ranges, 
                           const NanoFMI::Range &nr, bool split_all) {

    if (!nr.is_valid()) {
        return 0;
    }

    auto lb = next_ranges.lower_bound(nr);


    if (lb == next_ranges.end()) {
        next_ranges.insert(nr);
        return 1;
    }

    if (lb->same_range(nr)) {

        //PROBABILY THE PROBLEM
        //or maybe not
        if (!split_all && (nr.total_len_ > lb->total_len_ 
                            || (nr.prob_ > lb->prob_ && nr.total_len_ == lb->total_len_))) {
            Range old_q = *lb;

            //std::cout << nr.match_len_ << " vs " << lb->match_len_ << "\n";
            
            lb->parent_->child_count_--; //TODO: go further? 
            next_ranges.erase(lb); 
            next_ranges.insert(nr); //TODO: store lb next iter to speed up
            return 1;
        }
        if (nr.parent_ != NULL)
            nr.parent_->child_count_--;
        return 0; 
    }

    if (!split_all) {
        next_ranges.insert(nr);
        return 1;
    }

    auto start = lb;
    while (start != next_ranges.begin() && std::prev(start)->intersects(nr)){
        start--;
    }

    auto end = lb;
    while (end != next_ranges.end() && end->intersects(nr)) {
        end++;
    }

    if (start == next_ranges.end() || (start == lb && end == lb)) {
        next_ranges.insert(nr);
        return 1;
    }

    std::list<Range> split_ranges;

    Range rr(nr); //Right range

    for (auto pr = start; pr != end; pr++) {
        Range lr = rr.split_range(*pr); //Left range

        if (lr.is_valid()) 
            split_ranges.push_back(lr);

        if (!lr.is_valid())
            break;
    }

    if (rr.is_valid())
        split_ranges.push_back(rr);

    next_ranges.insert(split_ranges.begin(), split_ranges.end());

    return split_ranges.size();
}

bool NanoFMI::Range::intersects(const Range &q) const {
    return (start_ >= q.start_ && start_ <= q.end_) ||
           (end_ >= q.start_ && end_ <= q.end_);
}

NanoFMI::Range NanoFMI::Range::split_range(const NanoFMI::Range &r) { 

    //if (prob_ < q.prob_)
    //    return Range(fmi_);

    Range left(fmi_);
    if (start_ < r.start_) {
        left = Range(*this);
        left.end_ = r.start_ - 1;
    }

    if (end_ > r.end_) {
        start_ = r.end_ + 1;
    } else {
        start_ = end_ = prob_ = 0;
    }

    if (parent_ != NULL) {
        if (left.is_valid())
            left.parent_->child_count_++;

        if (!is_valid())
            parent_->child_count_--;
    }

    return left;
}


NanoFMI::Range& NanoFMI::Range::operator=(const NanoFMI::Range& prev) {
    fmi_ = prev.fmi_;
    k_id_ = prev.k_id_;
    event_ = prev.event_;
    start_ = prev.start_;
    end_ = prev.end_;
    prob_ = prev.prob_;
    match_len_ = prev.match_len_;
    stay_ = prev.stay_;
    parent_ = prev.parent_;
    //children_ = std::list<Range const *>(prev.children_);
    child_count_ = prev.child_count_;
    return *this;
}

NanoFMI::Range::Range(NanoFMI &fmi)
    : fmi_(fmi), 
      k_id_(0),
      event_(-1),
      start_(0), 
      end_(0), 
      match_len_(0),
      total_len_(0),
      prob_(0),
      stay_(false),
      parent_(NULL),
      child_count_(0) {}

//Initial match constructor
NanoFMI::Range::Range(NanoFMI &fmi, mer_id k_id, int event, double prob)
    : fmi_(fmi), 
      k_id_(k_id),
      event_(event),
      start_(fmi.mer_f_starts_[k_id]), 
      end_(fmi.mer_f_starts_[k_id] + fmi.mer_counts_[k_id] - 1), 
      match_len_(1),
      total_len_(1),
      prob_(prob),
      stay_(false),
      parent_(NULL),
      child_count_(0) {}

//"next" constructor
NanoFMI::Range::Range(const NanoFMI::Range &prev, mer_id k_id, double prob)
    : fmi_(prev.fmi_),
      k_id_(k_id),
      stay_(false),
      parent_(&prev) {

    int min = fmi_.get_tally(k_id, prev.start_ - 1),
        max = fmi_.get_tally(k_id, prev.end_);

    //Candidate k-mer occurs in the range
    if (min < max) {
        event_ = prev.event_ - 1,
        start_ = fmi_.mer_f_starts_[k_id] + min;
        end_ = fmi_.mer_f_starts_[k_id] + max - 1;
        match_len_ = prev.match_len_ + 1;
        total_len_ = prev.total_len_ + 1;
        prob_ = prev.prob_ + prob;
        prev.child_count_++;
        //prev.children_.push_back(this);
    } else {
        start_ = end_ = 0;

    }
}

//"stay" constructor
NanoFMI::Range::Range(const NanoFMI::Range &prev, double prob)
    : fmi_(prev.fmi_), 
      k_id_(prev.k_id_),
      event_(prev.event_ - 1),
      start_(prev.start_), 
      end_(prev.end_), 
      match_len_(prev.match_len_),
      total_len_(prev.total_len_+1),
      prob_(prev.prob_ + prob),
      stay_(true),
      parent_(&prev) {

        //prev.children_.push_back(this);
        prev.child_count_++;
          
    }

//Copy constructor
NanoFMI::Range::Range(const NanoFMI::Range &prev)
    : fmi_(prev.fmi_), 
      k_id_(prev.k_id_),
      event_(prev.event_),
      start_(prev.start_), 
      end_(prev.end_), 
      match_len_(prev.match_len_),
      total_len_(prev.total_len_),
      prob_(prev.prob_),
      stay_(prev.stay_),
      parent_(prev.parent_),
      child_count_(prev.child_count_) {
      //children_(std::list<Range const *>(prev.children_)) {
} 

NanoFMI::Range::~Range() {
    //delete parents_;
    //delete children_;
}

bool NanoFMI::Range::same_range(const NanoFMI::Range &q) const {
    return start_ == q.start_ && end_ == q.end_;
}

bool NanoFMI::Range::is_valid() const {
    return start_ <= end_ && (start_ != 0 || end_ != 0);
}


void NanoFMI::Range::print_info() const {
}


bool operator< (const NanoFMI::Range &q1, const NanoFMI::Range &q2) {
    if (q1.event_ > q2.event_) {
        return true;
        
    } else if (q1.event_ == q2.event_) {

        if (q1.start_ < q2.start_)
            return true;

        if (q1.start_ == q2.start_ && q1.end_ < q2.end_)
            return true;
    }
    
    return false;
}
