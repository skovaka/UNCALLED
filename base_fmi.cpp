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


BaseFMI::BaseFMI() {
    loaded_ = false;

}

//Reads a model directly from a file and creates the FM index from the given reference
BaseFMI::BaseFMI(std::string seq, unsigned int tally_dist) {

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

    //std::cerr << "SA init time: " << timer.lap() << "\n";

    //Create the suffix array
    std::sort(suffix_ar.begin(), suffix_ar.end(), (BaseFMI) *this);
    suffix_ar_.swap(suffix_ar);

    //std::cerr << "SA sort time: " << timer.lap() << "\n";

    //Allocate space for other data structures
    bwt_.resize(seq.size());
    f_starts_.resize(ALPH_SIZE);
    counts_.resize(ALPH_SIZE);
    tally_.resize(ALPH_SIZE);

    for (size_t i = 0; i < ALPH_SIZE; i++)
        tally_[i].resize((seq.size() / tally_dist_) + 1, -1);
    
    //std::cerr << "FM init time: " << timer.lap() << "\n";

    unsigned int tally_mod = tally_dist_;

    //Single pass to generate BWT and other datastructures
    for (size_t i = 0; i < suffix_ar_.size(); i++) {
        
        //Fill in BWT
        if (suffix_ar_[i] > 0) {
            bwt_[i] = char_to_base(seq[suffix_ar_[i]-1]);
        } else {
            bwt_[i] = char_to_base(seq.back());
        }

        //Update 6-mer counts
        if (bwt_[i] < ALPH_SIZE)
            counts_[bwt_[i]]++;
        
        //Update tally array
        if (tally_mod == tally_dist_) {
            for (size_t j = 0; j < ALPH_SIZE; j++) {
                tally_[j][i / tally_dist_] = counts_[j];
            }
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    //std::cerr << "FM build time: " << timer.lap() << "\n";
    
    //TODO: store as range?
    f_starts_[0] = 1;
    for (size_t i = 1; i < ALPH_SIZE; i++) {
        f_starts_[i] = f_starts_[i-1] + counts_[i-1];
    }

    
    //Fill in last entry in tally array if needed
    if (seq.size() % tally_dist_ == 0) {
        for (size_t i = 0; i < ALPH_SIZE; i++) {
            tally_[i][tally_[i].size()-1] = counts_[i];
        }
    }
    loaded_ = true;
}

BaseFMI::BaseFMI(std::ifstream &infile, unsigned int tally_dist) {
    tally_dist_ = tally_dist;
    
    size_t size = 0;
    
    infile >> size;

    std::cerr << size << "\n";
    
    suffix_ar_.resize(size);
    bwt_.resize(size);
    f_starts_.resize(ALPH_SIZE);
    counts_.resize(ALPH_SIZE);
    tally_.resize(ALPH_SIZE);

    for (size_t i = 0; i < ALPH_SIZE; i++)
        tally_[i].resize((size / tally_dist_) + 1, -1);
    
    unsigned int tally_mod = tally_dist_;
    unsigned int bwt_i;

    for (size_t i = 0; i < size; i++) {
        infile >> bwt_i >> suffix_ar_[i];

        bwt_[i] = (base_t) bwt_i;

        if (bwt_[i] < ALPH_SIZE)
            counts_[bwt_[i]]++;
        
        //Update tally array
        if (tally_mod == tally_dist_) {
            for (size_t j = 0; j < ALPH_SIZE; j++) {
                tally_[j][i / tally_dist_] = counts_[j];
            }
            tally_mod = 0;
        }
        tally_mod += 1;

        //if (++i_mod == update_interval) {
        //    std::cerr << (100.0 * i) / size << "%\n";
        //    i_mod = 0;
        //}
    }

    //TODO: store as range?
    f_starts_[0] = 1;
    for (size_t i = 1; i < ALPH_SIZE; i++) {
        f_starts_[i] = f_starts_[i-1] + counts_[i-1];
    }
    
    //Fill in last entry in tally array if needed
    if (size % tally_dist_ == 0) {
        for (size_t i = 0; i < ALPH_SIZE; i++) {
            tally_[i][tally_[i].size()-1] = counts_[i];
        }
    }
    loaded_ = true;
}

void BaseFMI::save(const std::string &filename) {
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

    base_t c1, c2;
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

Range BaseFMI::get_neighbor(Range range, base_t base) const {
    unsigned int min = get_tally(base, range.start_ - 1);
    unsigned int max = get_tally(base, range.end_);

    if (min >= 0 && max >= 0 && min < max) {
        unsigned int base_st = f_starts_[base];
        return Range(base_st + min, base_st + max - 1);
    }

    return Range();
}

unsigned int BaseFMI::get_tally(base_t base, unsigned int loc) const {
    if (loc < 0)
        return -1;

    //Closest checkpoint < i
    unsigned int cp = (loc / tally_dist_) * tally_dist_; 

    //Check if checkpoint after i is closer
    if (loc - cp > (cp + tally_dist_) - loc 
        && cp + tally_dist_ < bwt_.size()) {
        cp += tally_dist_;
    }

    base_t k = base;
    unsigned int count = tally_[k][cp / tally_dist_];

    //int cp_dist = cp - loc;

    if (cp > loc) {
        for (size_t i = loc + 1; i <= cp; i++) {
            count -= bwt_[i] == k;
        }

    } else if (cp < loc) {
        for (size_t i = loc; i > cp; i--) {
            count += bwt_[i] == k;
        }
    }
    
    return count;
}


//TODO: Maybe store f as ranges?
Range BaseFMI::get_full_range(base_t base) const {
    return Range(f_starts_[base], f_starts_[base] + counts_[base] - 1 );
}  

size_t BaseFMI::sa(size_t i) const {
    return suffix_ar_[i];
}

size_t BaseFMI::size() const {
    return bwt_.size();
}

//Range BaseFMI::get_kmer_range(const std::string &seq) const {
//    Range r = get_full_range(seq.back());
//
//    for (size_t i = seq.size()-2; i < seq.size(); i--) {
//        if (!r.is_valid()) {
//            break;
//        }
//
//        r = get_neighbor(r, seq[i]);
//    }
//
//    return r;
//}


