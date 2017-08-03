#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stddef.h>
#include <list>
#include <set>
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
    int df = pow(var1/e1.length + var2/e2.length, 2) / (pow(var1/e1.length, 2) / (e1.length-1) + pow(var2/e2.length, 2) / (e2.length-1));

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

bool operator< (const FMQuery &q1, const FMQuery &q2) {
    //return q1.start < q2.start || q1.end < q2.end || q1.avg_prob() < q2.avg_prob();
    return q1.start < q2.start ||
           q1.end < q2.end ||
           q1.avg_prob() > q2.avg_prob() ||
           q1.stays > q2.stays||
           q1.match_len < q2.match_len;
}

bool FMQuery::same_range(const FMQuery &q) const {
    return start == q.start && end == q.end && match_len == q.match_len;
}

float FMQuery::avg_prob() const {
    return prob_sum / (match_len + stays);
}

//Aligns the vector of events to the reference using LF mapping
//Returns the number of exact alignments
int NanoFMI::lf_map(std::vector<Event> &events, int seed_end, 
                    int seed_len, NormParams norm_params) {

    float EVENT_THRESH = 0.00025,
          SEED_THRESH  = 0.05,
          STAY_THRESH  = 0.01;

    //Stores ranges corrasponding to the F array and the L array (the BWT)
    std::vector< std::set<FMQuery> > f_locs(model_->kmer_count()), 
                                     l_locs(model_->kmer_count());
    std::list<mer_id> f_mers, l_mers;
    std::set<FMQuery> results;

    //Match the first event

    //f_locs = match_event(events.back(), scale);

    double prob;
    for (mer_id i = 0; i < f_locs.size(); i++) {
        prob = model_->event_match_prob(events[seed_end], i, norm_params);
        if (prob >= EVENT_THRESH) {
            f_locs[i].emplace(mer_f_starts_[i], mer_counts_[i], prob);
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
            //std::cout << "Fmer " << f << "\n";

            prob = model_->event_match_prob(events[i], f, norm_params);

            if (stay && prob >= EVENT_THRESH) {
                if (l_locs[f].empty())
                    l_mers.push_back(f);

                for (auto fq = f_locs[f].begin(); fq != f_locs[f].end(); ++fq) {
                    l_locs[f].emplace(*fq, prob);
                }
            }

            //Check neighboring k-mers
            //std::vector<mer_id> &neighbors = mer_neighbors[f];

            //TODO: switch this and innermost for loop?
            //compute tally for all neighbors at once, delete f_locs as we go
            for (auto n = model_->neighbor_begin(f); 
                    n != model_->neighbor_end(f); n++) {
                
                //TODO: Never compute same prob twice ??
                prob = model_->event_match_prob(events[i], *n, norm_params);

                if(prob >= EVENT_THRESH) {

                    int loop_count = 0;

                    for (auto fq = f_locs[f].begin(); fq != f_locs[f].end(); ++fq) {
                        int min_tally = get_tally(*n, (fq->start) - 1),
                            max_tally = get_tally(*n, fq->end);
                        
                        //Candidate k-mer occurs in the range
                        if (min_tally < max_tally) {
                            FMQuery lq(*fq, mer_f_starts_[*n]+min_tally, max_tally-min_tally, prob);

                            if ((fq->match_len)+1 < seed_len) {
                                if (l_locs[*n].empty())
                                    l_mers.push_back(*n);

                                l_locs[*n].insert(lq);

                                mer_matched = true;
                            } else {
                                results.insert(lq);
                            }
                        }
                        loop_count++;
                    }
                }
            }
            f_locs[f].clear();
        }

        //Update current ranges
        for (auto m = l_mers.begin(); m != l_mers.end(); ++m) {
            std::set<FMQuery>::iterator prev = f_locs[*m].end();
            while (!l_locs[*m].empty()) {
                std::set<FMQuery>::iterator q = l_locs[*m].begin();

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
    
    //Count up the number of alignments
    //Currently not doing anything with location of alignments
    int count = 0;
    
    if (!results.empty()) {
        std::set<FMQuery>::iterator prev = results.end();
        for (auto f_loc = results.begin(); f_loc != results.end(); ++f_loc) {
            if ((prev == results.end() || !f_loc->same_range(*prev))
                && f_loc->avg_prob() >= SEED_THRESH) {

                for (int s = f_loc->start; s <= f_loc->end; s++) {
                    std::cout << seed_end - f_loc->match_len - f_loc->stays + 1 << " " 
                              << suffix_ar_[s] << " " 
                              << f_loc->avg_prob() << " " 
                              << f_loc->match_len << " " 
                              << f_loc->stays << "\n";

                    //std::cout << seed_end - f_loc->match_len - f_loc->stays + 1 << " " << suffix_ar[s] << " " << f_loc->avg_prob() << "\n";
                        
                    //std::cout << "\n";
                    count += 1; 
                }
            }
            prev = f_loc;
        }
    }

    return count;
}



