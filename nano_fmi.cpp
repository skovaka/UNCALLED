#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stddef.h>
#include <list>
#include <set>
#include <utility>
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
                    std::vector<Event> &events, int seed_end, 
                    int seed_len, NormParams norm_params) {

    float EVENT_THRESH = 0.00025,
          SEED_THRESH  = 0.05,
          STAY_THRESH  = 0.01;

    //Stores ranges corrasponding to the F array and the L array (the BWT)
    std::vector< std::set<Query> > f_locs(model_->kmer_count()), 
                                   l_locs(model_->kmer_count());
    std::list<mer_id> f_mers, l_mers;
    std::set<Query> finished;

    //Match the first event

    double prob;
    for (mer_id i = 0; i < f_locs.size(); i++) {
        prob = model_->event_match_prob(events[seed_end], i, norm_params);
        if (prob >= EVENT_THRESH) {
            f_locs[i].emplace(*this, i, prob);
            f_mers.push_back(i);
        }
    }

    bool mer_matched = true,
         stay = false;

    //Match all other events backwards
    int i = seed_end - 1;
    while (i >= 0 && mer_matched) {
        mer_matched = false;

        //stay = get_stay_prob(events[i], events[i+1]) >= STAY_THRESH;
        stay = true;

        //Find all candidate k-mers
        while (!f_mers.empty()) {
            mer_id f = f_mers.front();
            f_mers.pop_front();

            prob = model_->event_match_prob(events[i], f, norm_params);

            if (stay && prob >= EVENT_THRESH) {
                if (l_locs[f].empty())
                    l_mers.push_back(f);

                for (auto fq = f_locs[f].begin(); fq != f_locs[f].end(); ++fq) {
                    l_locs[f].emplace(*fq, prob);
                }
            }

            auto neighbors = model_->get_neighbors(f);
            for (auto n = neighbors.first; n != neighbors.second; n++) {
                
                //TODO: Never compute same prob twice ??
                prob = model_->event_match_prob(events[i], *n, norm_params);

                if(prob >= EVENT_THRESH) {

                    for (auto fq = f_locs[f].begin(); fq != f_locs[f].end(); ++fq) {

                        Query lq(*fq, *n, prob);

                        if (lq.is_valid()) {
                            if (lq.match_len_ < seed_len) {
                                if (l_locs[*n].empty())
                                    l_mers.push_back(*n);

                                l_locs[*n].insert(lq);

                                mer_matched = true;
                            } else {
                                finished.insert(lq);
                            }
                        }
                    }
                }
            }
            f_locs[f].clear();
        }

        //Update current ranges
        for (auto m = l_mers.begin(); m != l_mers.end(); ++m) {
            std::set<Query>::iterator prev = f_locs[*m].end();
            while (!l_locs[*m].empty()) {
                std::set<Query>::iterator q = l_locs[*m].begin();

                if (prev == f_locs[*m].end() || !q->same_range(*prev))
                    prev = f_locs[*m].insert(prev, *q);
    
                l_locs[*m].erase(q);
            }
        }
        l_mers.swap(f_mers);
        
        //Stop search if not matches found
        if (!mer_matched)
            break;

        i--;
    }
    

    std::vector<Result> results;

    auto prev = finished.end();
    for (auto qry = finished.begin(); qry != finished.end(); qry++) {
        if ( prev == finished.end() || !qry->same_range(*prev) ) 
            qry->add_results(results, seed_end, SEED_THRESH);
        prev = qry;
    }

    return results;
}

bool operator< (const NanoFMI::Query &q1, const NanoFMI::Query &q2) {
    //return q1.start < q2.start || q1.end < q2.end || q1.avg_prob() < q2.avg_prob();
    return q1.start_ < q2.start_ ||
           q1.end_ < q2.end_ ||
           q1.avg_prob() > q2.avg_prob() ||
           q1.stays_ > q2.stays_ ||
           q1.match_len_ < q2.match_len_;
}

//Initial match constructor
NanoFMI::Query::Query(NanoFMI &fmi, mer_id k_id, float prob)
    : fmi_(fmi), 
      start_(fmi.mer_f_starts_[k_id]), 
      end_(fmi.mer_f_starts_[k_id] + fmi.mer_counts_[k_id] - 1), 
      match_len_(1), 
      stays_(0), 
      prob_sum_(prob) {}

//"next" constructor
NanoFMI::Query::Query(const NanoFMI::Query &prev, mer_id k_id, double prob)
    : fmi_(prev.fmi_) {
    int min = fmi_.get_tally(k_id, prev.start_ - 1),
        max = fmi_.get_tally(k_id, prev.end_);
    
    //Candidate k-mer occurs in the range
    if (min < max) {
        start_ = fmi_.mer_f_starts_[k_id] + min;
        end_ = fmi_.mer_f_starts_[k_id] + max - 1;
        match_len_ = prev.match_len_ + 1;
        stays_ = prev.stays_;
        prob_sum_ = prev.prob_sum_ + prob;
    } else {
        start_ = end_ = match_len_ = stays_ = prob_sum_ = 0;
    }
}

//"stay" constructor
NanoFMI::Query::Query(const NanoFMI::Query &prev, float prob)
    : fmi_(prev.fmi_), 
      start_(prev.start_), 
      end_(prev.end_), 
      match_len_(prev.match_len_), 
      stays_(prev.stays_ + 1), 
      prob_sum_(prev.prob_sum_ + prob) {}

//Copy constructor
NanoFMI::Query::Query(const NanoFMI::Query &prev)
    : fmi_(prev.fmi_), 
      start_(prev.start_), 
      end_(prev.end_), 
      match_len_(prev.match_len_), 
      stays_(prev.stays_), 
      prob_sum_(prev.prob_sum_) {} 

bool NanoFMI::Query::add_results(std::vector<NanoFMI::Result> &results, 
                                 int query_end, double min_prob) const {
    if (avg_prob() >= min_prob) {
        Result r;
        for (int s = start_; s <= end_; s++) {
            r.qry_start = query_end - match_len_ - stays_ + 1;
            r.qry_end = query_end;
            r.ref_start = fmi_.suffix_ar_[s];
            r.ref_end = fmi_.suffix_ar_[s] + match_len_ - 1;
            r.prob = avg_prob();
            results.push_back(r);
        }
        return true;
    }
    return false;
}

bool NanoFMI::Query::same_range(const NanoFMI::Query &q) const {
    return start_ == q.start_ && end_ == q.end_ && match_len_ == q.match_len_;
}

bool NanoFMI::Query::is_valid() {
    return start_ != 0 || end_ != 0 || match_len_ != 0 || 
           stays_ != 0 || prob_sum_ != 0;
}

float NanoFMI::Query::avg_prob() const {
    return prob_sum_ / (match_len_ + stays_);
}


