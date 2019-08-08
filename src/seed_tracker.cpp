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
#include "params.hpp"

#define DEBUG

SeedGroup::SeedGroup() 
    : evt_st_(1),
      evt_en_(0),
      total_len_(0) {
}

SeedGroup::SeedGroup(Range ref_st, u32 evt_st)
    : ref_st_(ref_st.start_),
      ref_en_(ref_st),
      evt_st_(evt_st),
      evt_en_(evt_st),
      total_len_(ref_st.length()) {}

SeedGroup::SeedGroup(const SeedGroup &r)
    : ref_st_(r.ref_st_),
      ref_en_(r.ref_en_),
      evt_st_(r.evt_st_),
      evt_en_(r.evt_en_),
      total_len_(r.total_len_) {}


u8 SeedGroup::update(SeedGroup &new_aln) {
    u8 growth = 0;
    if (new_aln.ref_en_.start_ < ref_en_.end_) {
        if (new_aln.ref_en_.end_ > ref_en_.end_) {
            growth = new_aln.ref_en_.end_ - ref_en_.end_;
            ref_en_ = new_aln.ref_en_;
        } else {
            ref_en_.start_ = new_aln.ref_en_.start_;
        }
    } else {
        growth = new_aln.total_len_;
        ref_en_ = new_aln.ref_en_;
    }

    evt_en_ = new_aln.evt_en_;
    total_len_ += growth;
    return growth;
}

Range SeedGroup::ref_range() const {
    return Range(ref_st_, ref_en_.end_);
}

void SeedGroup::print(std::ostream &out, bool newline = false, bool print_all = false) const {
    out << total_len_ << "\t";

    out << ref_st_;

    out << "-" << ref_en_.end_ << "\t" 
               << evt_st_ << "-" 
               << evt_en_;

    if (newline)
        out << "\n";
}

bool SeedGroup::is_valid() {
    return evt_st_ <= evt_en_;
}


bool operator< (const SeedGroup &r1, const SeedGroup &r2) {
    if (r1.ref_en_.start_ != r2.ref_en_.start_)
        return r1.ref_en_.start_ > r2.ref_en_.start_;

    return r1.evt_en_ > r2.evt_en_;
}

std::ostream &operator<< (std::ostream &out, const SeedGroup &a) {
    out << a.ref_st_ << "-" << a.ref_en_.end_ << "\t"
        << a.evt_st_ << "-" << (a.evt_en_) << "\t"
        << a.total_len_;
    return out;
}

SeedTracker::SeedTracker() {
    reset();
}

void SeedTracker::reset() {
    alignments_.clear();
    all_lens_.clear();
    max_map_ = NULL_ALN;
    len_sum_ = 0;
}

SeedGroup SeedTracker::get_final() {
    if (max_map_.total_len_ < PARAMS.min_aln_len || 
        all_lens_.size() < 2) return NULL_ALN;

    float mean_len = len_sum_ / alignments_.size();
    float second_len = *std::next(all_lens_.rbegin());

    if (PARAMS.check_map_conf(max_map_.total_len_, mean_len, second_len)) {

        //print(std::cout, 10);
        return max_map_;
    }
    
    return NULL_ALN;
}

SeedGroup SeedTracker::get_best() {
    return max_map_;
}

float SeedTracker::get_top_conf() {
    return (float) max_map_.total_len_ / (*std::next(all_lens_.rbegin()));
}

float SeedTracker::get_mean_conf() {
    return max_map_.total_len_ / (len_sum_ / alignments_.size());
}

void SeedTracker::add_seed(u64 ref_en, u32 ref_len, u32 evt_st) {

    SeedGroup new_aln(Range(ref_en-ref_len+1, ref_en), evt_st);
    //new_aln.print(std::cout, true, false);

    //Locations sorted by decreasing ref_en_.start
    //Find the largest aln s.t. aln->ref_en_.start <= new_aln.ref_en_.start
    //AKA r1 <= r2
    auto aln = alignments_.lower_bound(new_aln),
         aln_match = alignments_.end();

    u64 e2 = new_aln.evt_en_, //new event aln
        r2 = new_aln.ref_en_.start_; //new ref aln

    while (aln != alignments_.end()) {
        u64 e1 = aln->evt_en_, //old event aln
            r1 = aln->ref_en_.start_; //old ref aln

        //We know r1 <= r2 because of alnation sort order

        bool higher_sup = aln_match == alignments_.end() 
                       || aln_match->total_len_ < aln->total_len_,
             
             in_range = e1 <= e2 && //event aln must increase
                        r2 - r1 <= e2 - e1 && //evt increases more than ref (+ skip)
                        (r2 - r1) >= (e2 - e1) / 12; //evt doesn't increase too much
             
        if (higher_sup && in_range) {
            aln_match = aln;
        } else if (r2 - r1 >= e2) {
            break;
        }

        aln++;
    }

    if (aln_match != alignments_.end()) {
        SeedGroup a = *aln_match;

        u32 prev_len = a.total_len_;
        a.update(new_aln);

        if (a.total_len_ != prev_len) {
            len_sum_ += a.total_len_ - prev_len;
            auto l = all_lens_.find(prev_len);
            all_lens_.insert(l, a.total_len_);
            all_lens_.erase(l);

            if (a.total_len_ >= PARAMS.min_aln_len && a.total_len_ > max_map_.total_len_) {
                max_map_ = a;
            }
        }

        auto hint = std::next(aln_match);
        alignments_.erase(aln_match);
        alignments_.insert(hint, a);

    }

    alignments_.insert(new_aln);
    all_lens_.insert(new_aln.total_len_);
    len_sum_ += new_aln.total_len_;

    if (new_aln.total_len_ >= PARAMS.min_aln_len && new_aln.total_len_ > max_map_.total_len_) {
        max_map_ = new_aln;
    }
}

void SeedTracker::print(std::ostream &out, u16 max_out = 10) {
    if (alignments_.empty()) {
        return;
    }

    std::vector<SeedGroup> alns_sort(alignments_.begin(),
                                   alignments_.end());

    std::sort(alns_sort.begin(), alns_sort.end(),
              [](const SeedGroup &a, const SeedGroup &b) -> bool {
                  return a.total_len_ > b.total_len_;
              });

    Range top_ref = alns_sort[0].ref_range();
    float top_len = alns_sort[0].total_len_;

    for (unsigned int i = 0; i < min(max_out, alns_sort.size()); i++) {
        float overlap = top_ref.get_recp_overlap(alns_sort[i].ref_range()),
               len_ratio = top_len / alns_sort[i].total_len_;

        alns_sort[i].print(out, false);
        out << "\t" << len_ratio << "\t" << overlap << "\n";
    }
}

