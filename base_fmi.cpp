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
#include "base_fmi.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)

int max(int a, int b) {
    return a > b ? a : b;
}

int min(int a, int b) {
    return a < b ? a : b;
}

char base_to_idx(char base) {
    switch (base) {
        case 'A':
        return 0;
        case 'C':
        return 1;
        case 'G':
        return 2;
        case 'T':
        return 3;
        default:
        return 4;
    }
}

BaseFMI::BaseFMI() {
    loaded_ = false;

}

//Reads a model directly from a file and creates the FM index from the given reference
BaseFMI::BaseFMI(std::string seq, int tally_dist) {

    seq_ = &seq;
    tally_dist_ = tally_dist;

    //For outputting time
    Timer timer;

    //Init suffix array
    //Not using suffix_ar instance var speeds up sorting significantly
    std::vector<unsigned int> suffix_ar(seq.size());
    for (unsigned int i = 0; i < suffix_ar.size(); i++) {
        suffix_ar[i] = i;
    }

    std::cerr << "SA init time: " << timer.lap() << "\n";

    //Create the suffix array
    std::sort(suffix_ar.begin(), suffix_ar.end(), *this);
    suffix_ar_.swap(suffix_ar);

    std::cerr << "SA sort time: " << timer.lap() << "\n";

    //Allocate space for other data structures
    bwt_ = std::string(seq.size(), '$');
    f_starts_.resize(4);
    counts_.resize(4);
    tally_.resize(4);

    count_tmp_.resize(4, 0);

    for (size_t i = 0; i < 4; i++)
        tally_[i].resize((seq.size() / tally_dist_) + 1, -1);
    
    std::cerr << "FM init time: " << timer.lap() << "\n";

    int tally_mod = tally_dist_;

    //Single pass to generate BWT and other datastructures
    for (unsigned int i = 0; i < suffix_ar_.size(); i++) {
        
        //Fill in BWT
        if (suffix_ar_[i] > 0) {
            bwt_[i] = base_to_idx(seq[suffix_ar_[i]-1]);
        } else {
            bwt_[i] = base_to_idx(seq.back());
        }


        //Update 6-mer counts
        counts_[bwt_[i]]++;
        
        //Update tally array
        if (tally_mod == tally_dist_) {
            for (size_t j = 0; j < 4; j++) {
                tally_[j][i / tally_dist_] = counts_[j];
            }
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    std::cerr << "FM build time: " << timer.lap() << "\n";
    
    //TODO: store as range?
    f_starts_[0] = 1;
    for (size_t i = 1; i < 4; i++) {
        f_starts_[i] = f_starts_[i-1] + counts_[i-1];
    }

    
    //Fill in last entry in tally array if needed
    if (seq.size() % tally_dist_ == 0) {
        for (size_t i = 0; i < 4; i++) {
            tally_[i][tally_[i].size()-1] = counts_[i];
        }
    }
    loaded_ = true;
}

BaseFMI::BaseFMI(std::ifstream &infile, int tally_dist) {
    tally_dist_ = tally_dist;
    
    size_t size = 0;
    
    infile >> size;

    std::cerr << size << "\n";
    
    suffix_ar_.resize(size);
    bwt_ = std::string(size, '$');
    f_starts_.resize(4);
    counts_.resize(4);
    tally_.resize(4);
    count_tmp_.resize(4, 0);

    for (size_t i = 0; i < 4; i++)
        tally_[i].resize((size / tally_dist_) + 1, -1);
    
    int tally_mod = tally_dist_;
    int bwt_i;

    for (size_t i = 0; i < size; i++) {
        infile >> bwt_i >> suffix_ar_[i];

        bwt_[i] = (char) bwt_i;

        if (bwt_[i] < 4)
            counts_[bwt_[i]]++;
        
        //Update tally array
        if (tally_mod == tally_dist_) {
            for (size_t j = 0; j < 4; j++) {
                tally_[j][i / tally_dist_] = counts_[j];
            }
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    //TODO: store as range?
    f_starts_[0] = 1;
    for (size_t i = 1; i < 4; i++) {
        f_starts_[i] = f_starts_[i-1] + counts_[i-1];
    }
    
    //Fill in last entry in tally array if needed
    if (size % tally_dist_ == 0) {
        for (size_t i = 0; i < 4; i++) {
            tally_[i][tally_[i].size()-1] = counts_[i];
        }
    }
    loaded_ = true;
}

void BaseFMI::save(std::string filename) {
    std::ofstream out(filename);

    out << bwt_.size() << "\n";
    for (size_t i = 0; i < bwt_.size(); i++) {
        out << (int) bwt_[i] << "\t" << suffix_ar_[i] << "\n";
    }

    out.close();
}

//Returns true if the suffix of *mer_seq_tmp starting at rot1 is less than that
//starting at rot2. Used to build suffix array.
bool BaseFMI::operator() (unsigned int rot1, unsigned int rot2) {

    int c1, c2;
    for (unsigned int i = 0; i < seq_->size(); i++) {
        
        c1 = seq_->at(rot1 + i);
        c2 = seq_->at(rot2 + i);

        if (c1 == c2)
            continue;

        if (c2 == '$')
            return false;
        
        if (c1 == '$' || c1 < c2)
            return true;

       return false;
    }

    return false;
}

Range BaseFMI::get_neighbor(Range range, char base) const {
    int min = get_tally(base, range.start_ - 1);
    int max = get_tally(base, range.end_);

    if (min >= 0 && max >= 0 && min < max) {
        int base_st = f_starts_[base_to_idx(base)];
        return Range(base_st + min, base_st + max - 1);
    }

    return Range();
}

int BaseFMI::get_tally(char base, int loc) const {
    if (loc < 0)
        return -1;

    //Closest checkpoint < i
    int cp = (loc / tally_dist_) * tally_dist_; 

    //Check if checkpoint after i is closer
    if (loc - cp > (cp + tally_dist_) - loc 
            && cp + (unsigned) tally_dist_ < bwt_.size())
        cp += tally_dist_;

    int cp_dist = cp - loc; //TODO: just use cp

    int k = base_to_idx(base);

    int count = tally_[k][(loc + cp_dist) / tally_dist_];

    if (cp_dist > 0) {
        for (int i = loc+1; i <= loc + cp_dist; i++) {
            count -= bwt_[i] == k;
        }
    } else if (cp_dist < 0) {
        for (int i = loc; i > loc + cp_dist; i--) {
            count += bwt_[i] == k;
        }
    }
    
    return count;
}

//Returns the number of occurences of the given k-mer in the BWT up to and
//including the given index
std::list<int> BaseFMI::get_tallies(std::list<char> bases, int loc) const {
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

    for (auto k = bases.begin(); k != bases.end(); k++)
        count_tmp_[base_to_idx(*k)] = tally_[base_to_idx(*k)][(loc + cp_dist) / tally_dist_];

    if (cp_dist > 0) {
        for (int i = loc+1; i <= loc + cp_dist; i++) {
            count_tmp_[bwt_[i]]--;
        }
    } else if (cp_dist < 0) {
        for (int i = loc; i > loc + cp_dist; i--) {
            count_tmp_[bwt_[i]]++;
        }
    }
    
    for (auto k = bases.begin(); k != bases.end(); k++) {
        tallies.push_back(count_tmp_[base_to_idx(*k)]);
    }

    return tallies;
}


std::list<Range> BaseFMI::get_neigbhors(Range range, std::list<char> bases) const {
    std::list<Range> results;

    std::list<int> mins = get_tallies(bases, range.start_ - 1);
    std::list<int> maxs = get_tallies(bases, range.end_);

    auto base = bases.begin();
    auto min = mins.begin();
    auto max = maxs.begin();
        
    while (base != bases.end()) {
        if (*min < *max) {
            int base_st = f_starts_[base_to_idx(*base)];
            results.push_back(Range(base_st + *min, base_st + *max - 1));
        } else {
            results.push_back(Range());
        }

        base++;
        min++;
        max++;
    }

    return results;
}

//TODO: Maybe store f as ranges?
Range BaseFMI::get_full_range(char base) const {
    return Range(f_starts_[base_to_idx(base)], f_starts_[base_to_idx(base)] + counts_[base_to_idx(base)] -1 );
}

Range BaseFMI::get_kmer_range(const std::string &seq) const {
    Range r = get_full_range(seq.back());

    for (size_t i = seq.size()-2; i < seq.size(); i--) {
        if (!r.is_valid()) {
            break;
        }

        r = get_neighbor(r, seq[i]);
    }

    return r;
}

Range::Range(const Range &prev)
    : start_(prev.start_), 
      end_(prev.end_) {}

Range::Range(int start, int end) : start_(start), end_(end) {}

Range::Range() : start_(1), end_(0) {}

bool Range::intersects(const Range &q) const {
    return !(start_ > q.end_ || end_ < q.start_) &&
           !(q.start_ > end_ || q.end_ < start_);
}

int Range::length() const {
    return end_ - start_ + 1;
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

Range Range::intersect(const Range &r) const {
    if (!intersects(r)) {
        return Range();
    }
    
    return Range(max(start_, r.start_), min(end_, r.end_));
}

Range Range::merge(const Range &r) const {
    if (!intersects(r)) {
        return Range();
    }

    return Range(min(start_, r.start_), max(end_, r.end_));
}

double Range::get_recp_overlap(const Range &r) const {
    if (!intersects(r)) {
        return 0;
    }

    return double(intersect(r).length()) / double(merge(r).length());
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



