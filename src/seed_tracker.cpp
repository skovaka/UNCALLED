/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <set>
#include "seed_tracker.hpp"

const SeedTracker::Params SeedTracker::PRMS_DEF = {
    min_map_len   : 25,
    min_mean_conf : 6.00,
    min_top_conf  : 1.85
};

SeedCluster::SeedCluster() 
    : evt_st_(1),
      evt_en_(0),
      total_len_(0) {
}

SeedCluster::SeedCluster(Range ref_st, u32 evt_st)
    : ref_st_(ref_st.start_),
      ref_en_(ref_st),
      evt_st_(evt_st),
      evt_en_(evt_st),
      total_len_(ref_st.length()) {
}

//SeedCluster::SeedCluster(const SeedCluster &r)
//    : ref_st_(r.ref_st_),
//      ref_en_(r.ref_en_),
//      evt_st_(r.evt_st_),
//      evt_en_(r.evt_en_),
//      total_len_(r.total_len_) {}
//

u8 SeedCluster::update(SeedCluster &new_seed) {
    u8 growth = 0;
    if (new_seed.ref_en_.start_ < ref_en_.end_) {
        if (new_seed.ref_en_.end_ > ref_en_.end_) {
            growth = new_seed.ref_en_.end_ - ref_en_.end_;
            ref_en_ = new_seed.ref_en_;
        } else {
            ref_en_.start_ = new_seed.ref_en_.start_;
        }
    } else {
        growth = new_seed.total_len_;
        ref_en_ = new_seed.ref_en_;
    }

    evt_en_ = new_seed.evt_en_;
    total_len_ += growth;
    return growth;
}

Range SeedCluster::ref_range() const {
    return Range(ref_st_, ref_en_.end_);
}

void SeedCluster::print(std::ostream &out, bool newline = false, bool print_all = false) const {
    out << total_len_ << "\t";

    out << ref_st_;

    out << "-" << ref_en_.end_ << "\t" 
               << evt_st_ << "-" 
               << evt_en_;

    if (newline)
        out << "\n";
}

bool SeedCluster::is_valid() {
    return evt_st_ <= evt_en_;
}


bool operator< (const SeedCluster &r1, const SeedCluster &r2) {
    if (r1.ref_en_.start_ != r2.ref_en_.start_)
        return r1.ref_en_.start_ > r2.ref_en_.start_;

    return r1.evt_en_ > r2.evt_en_;
}

std::ostream &operator<< (std::ostream &out, const SeedCluster &a) {
    out << a.ref_st_ << "-" << a.ref_en_.end_ << "\t"
        << a.evt_st_ << "-" << (a.evt_en_) << "\t"
        << a.total_len_;
    return out;
}

SeedTracker::SeedTracker() : SeedTracker(PRMS_DEF) {}

SeedTracker::SeedTracker(Params prms) :
    PRMS(prms) {
    reset();
}

void SeedTracker::reset() {
    seed_clusters_.clear();
    all_lens_.clear();
    max_map_ = NULL_ALN;
    len_sum_ = 0;
}

bool SeedTracker::empty() {
    return seed_clusters_.empty();
}

SeedCluster SeedTracker::get_final() {
    if (max_map_.total_len_ < PRMS.min_map_len || 
        all_lens_.size() < 2) return NULL_ALN;

    float mean_len = len_sum_ / seed_clusters_.size();
    float second_len = *std::next(all_lens_.rbegin());

    if (check_map_conf(max_map_.total_len_, mean_len, second_len)) {

        //print(std::cout, 10);
        return max_map_;
    }
    
    return NULL_ALN;
}

SeedCluster SeedTracker::get_best() {
    return max_map_;
}

float SeedTracker::get_top_conf() {
    return (float) max_map_.total_len_ / (*std::next(all_lens_.rbegin()));
}

float SeedTracker::get_mean_conf() {
    return max_map_.total_len_ / (len_sum_ / seed_clusters_.size());
}

const SeedCluster &SeedTracker::add_seed(u64 ref_en, u32 ref_len, u32 evt_st) {
    SeedCluster new_seed(Range(ref_en-ref_len+1, ref_en), evt_st);
    
    //Locations sorted by decreasing ref_en_.start
    //Find the largest loc s.t. loc->ref_en_.start <= new_seed.ref_en_.start
    //AKA r1 <= r2
    auto loc = seed_clusters_.lower_bound(new_seed),
         loc_match = seed_clusters_.end();

    u64 e2 = new_seed.evt_en_, //new event loc
        r2 = new_seed.ref_en_.start_; //new ref loc

    while (loc != seed_clusters_.end()) {
        u64 e1 = loc->evt_en_, //old event loc
            r1 = loc->ref_en_.start_; //old ref loc

        //We know r1 <= r2 because of location sort order

        bool higher_sup = loc_match == seed_clusters_.end() 
                       || loc_match->total_len_ < loc->total_len_,
             
             in_range = e1 <= e2 && //event coord must increase
                        //r1 <= r2 &&
                        r2 - r1 <= e2 - e1 && //evt increases more than ref (+ skip)
                        
                        (r2 - r1) >= (e2 - e1) / 12; //evt doesn't increase too much
             
        if (higher_sup && in_range) {
            loc_match = loc;
        } else if (r2 - r1 >= e2) {
            break;
        }

        loc++;
    }

    auto ret = seed_clusters_.end();

    //If we find a matching seed cluster to join
    if (loc_match != seed_clusters_.end()) {
        SeedCluster a = *loc_match;

        u32 prev_len = a.total_len_;
        a.update(new_seed);

        if (a.total_len_ != prev_len) {
            len_sum_ += a.total_len_ - prev_len;
            auto l = all_lens_.find(prev_len);
            all_lens_.insert(l, a.total_len_);
            all_lens_.erase(l);

            if (a.total_len_ >= PRMS.min_map_len && a.total_len_ > max_map_.total_len_) {
                max_map_ = a;
            }
        }

        auto hint = std::next(loc_match);
        seed_clusters_.erase(loc_match);
        ret = seed_clusters_.insert(hint, a);
    } else {

        all_lens_.insert(new_seed.total_len_);
        len_sum_ += new_seed.total_len_;

        if (new_seed.total_len_ >= PRMS.min_map_len && new_seed.total_len_ > max_map_.total_len_) {
            max_map_ = new_seed;
        }

        #ifdef DEBUG_SEEDS
        new_seed.id_ = static_cast<u32>(seed_clusters_.size());
        #endif
        ret = seed_clusters_.insert(new_seed).first;
    }

    return *ret;
}

void SeedTracker::print(std::ostream &out, u16 max_out = 10) {
    if (seed_clusters_.empty()) {
        return;
    }

    std::vector<SeedCluster> seeds_sort(seed_clusters_.begin(),
                                     seed_clusters_.end());

    std::sort(seeds_sort.begin(), seeds_sort.end(),
              [](const SeedCluster &a, const SeedCluster &b) -> bool {
                  return a.total_len_ > b.total_len_;
              });

    Range top_ref = seeds_sort[0].ref_range();
    float top_len = seeds_sort[0].total_len_;

    for (unsigned int i = 0; i < std::min(max_out, (u16) seeds_sort.size()); i++) {
        float overlap = top_ref.get_recp_overlap(seeds_sort[i].ref_range()),
               len_ratio = top_len / seeds_sort[i].total_len_;

        seeds_sort[i].print(out, false);
        out << "\t" << len_ratio << "\t" << overlap << "\n";
    }
}

bool SeedTracker::check_map_conf(u32 seed_len, float mean_len, float second_len) {
    return (PRMS.min_mean_conf > 0 && seed_len / mean_len >= PRMS.min_mean_conf) ||
           (PRMS.min_top_conf > 0  && seed_len / second_len >= PRMS.min_top_conf);
}
