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

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)

//Reads a model directly from a file and creates the FM index from the given reference
NanoFMI::NanoFMI(int alph_size, std::vector<mer_id> &mer_seq, int tally_dist) {

    alph_size_ = alph_size;
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

    //Allocate space for other data structures
    bwt_.resize(mer_seq.size());
    mer_f_starts_.resize(alph_size_);
    mer_counts_.resize(alph_size_);
    mer_tally_.resize(alph_size_);

    mer_count_tmp_.resize(alph_size_, 0);

    for (mer_id i = 0; i < alph_size_; i++)
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
            for (mer_id j = 0; j < alph_size_; j++)
                mer_tally_[j][i / tally_dist_] = mer_counts_[j];
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    std::cerr << "FM build time: " << timer.lap() << "\n";
    
    //TODO: store as range?
    mer_f_starts_[0] = 1;
    for (mer_id i = 1; i < alph_size_; i++) {
        mer_f_starts_[i] = mer_f_starts_[i-1] + mer_counts_[i-1];
    }
    
    //Fill in last entry in tally array if needed
    if (mer_seq.size() % tally_dist_ == 0)
        for (mer_id i = 0; i < alph_size_; i++)
            mer_tally_[i][mer_tally_[i].size()-1] = mer_counts_[i];

}

//Returns true if the suffix of *mer_seq_tmp starting at rot1 is less than that
//starting at rot2. Used to build suffix array.
bool NanoFMI::operator() (unsigned int rot1, unsigned int rot2) {

    int c1, c2;
    for (unsigned int i = 0; i < mer_seq_->size(); i++) {
        
        c1 = mer_seq_->at(rot1 + i);
        c2 = mer_seq_->at(rot2 + i);

        if (c1 == c2)
            continue;

        if (c2 == alph_size_)
            return false;
        
        if (c1 == alph_size_ || c1 < c2)
            return true;

       return false;
    }

    return false;
}

//TODO: move to model
float NanoFMI::get_stay_prob(Event e1, Event e2) const {
    double var1 = e1.stdv*e1.stdv, var2 = e2.stdv*e2.stdv;

    double t = (e1.mean - e2.mean) / sqrt(var1/e1.length + var2/e2.length);

    int df = pow(var1/e1.length + var2/e2.length, 2) 
               / (pow(var1/e1.length, 2) / (e1.length-1) 
                    + pow(var2/e2.length, 2) / (e2.length-1));

    boost::math::students_t dist(df);
    double q = boost::math::cdf(boost::math::complement(dist, fabs(t)));

    return q;
}

//Returns the number of occurences of the given k-mer in the BWT up to and
//including the given index
std::list<int> NanoFMI::get_tallies(std::list<mer_id> kmers, int loc) const {
    std::list<int> tallies;
    if (loc < 0)
        return tallies;

    //Closest checkpoint < i
    int cp = (loc / tally_dist_) * tally_dist_; 

    //Check if checkpoint after i is closer
    if (loc - cp > (cp + tally_dist_) - loc 
            && cp + (unsigned) tally_dist_ < bwt_.size())
        cp += tally_dist_;

    int cp_dist = cp - loc; //TODO: just use cp

    for (auto k = kmers.begin(); k != kmers.end(); k++)
        mer_count_tmp_[*k] = mer_tally_[*k][(loc + cp_dist) / tally_dist_];

    if (cp_dist > 0) 
        for (int i = loc+1; i <= loc + cp_dist; i++)
            mer_count_tmp_[bwt_[i]]--;

    else if (cp_dist < 0)
        for (int i = loc; i > loc + cp_dist; i--)
            mer_count_tmp_[bwt_[i]]++;
    
    for (auto k = kmers.begin(); k != kmers.end(); k++) {
        tallies.push_back(mer_count_tmp_[*k]);
    }

    return tallies;
}

std::list<Range> NanoFMI::get_neigbhors(Range range, std::list<mer_id> kmers) const {
    std::list<Range> results;

    std::list<int> mins = get_tallies(kmers, range.start_ - 1);
    std::list<int> maxs = get_tallies(kmers, range.end_);

    auto kmer = kmers.begin();
    auto min = mins.begin();
    auto max = maxs.begin();
        
    while (kmer != kmers.end()) {
        if (*min < *max) {
            int kmer_st = mer_f_starts_[*kmer];
            results.push_back(Range(kmer_st + *min, kmer_st + *max - 1));
        } else {
            results.push_back(Range());
        }

        kmer++;
        min++;
        max++;
    }

    return results;
}

//TODO: Maybe store f as ranges?
Range NanoFMI::get_full_range(mer_id kmer) const {
    return Range(mer_f_starts_[kmer], mer_f_starts_[kmer] + mer_counts_[kmer] -1 );
}

Range::Range(const Range &prev)
    : start_(prev.start_), 
      end_(prev.end_) {}

Range::Range(int start, int end) : start_(start), end_(end) {}

Range::Range() : start_(1), end_(0) {}

bool Range::intersects(const Range &q) const {
    return (start_ >= q.start_ && start_ <= q.end_) ||
           (end_ >= q.start_ && end_ <= q.end_);
}

Range Range::split_range(const Range &r) { 

    Range left;
    if (start_ < r.start_) {
        left = Range(*this);
        left.end_ = r.start_ - 1;
    }

    if (end_ > r.end_) {
        start_ = r.end_ + 1;
    } else {
        start_ = 1;
        end_ = 0;
    }

    return left;
}

Range& Range::operator=(const Range& prev) {
    start_ = prev.start_;
    end_ = prev.end_;
    return *this;
}


bool Range::same_range(const Range &q) const {
    return start_ == q.start_ && end_ == q.end_;
}

bool Range::is_valid() const {
    return start_ <= end_;
}


bool operator< (const Range &q1, const Range &q2) {
    if (q1.start_ < q2.start_)
        return true;

    if (q1.start_ == q2.start_ && q1.end_ < q2.end_)
        return true;
    
    return false;
}



