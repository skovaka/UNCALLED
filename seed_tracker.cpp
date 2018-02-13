#include "seed_tracker.hpp"
#include "kmer_fmi.hpp"
#include <iostream>
#include <set>

//#define DEBUG

ReadAln::ReadAln() : evt_st_(-1) {
}

ReadAln::ReadAln(Range ref_en, int evt_en)
    : ref_st_(ref_en),
      evt_st_(evt_en),
      ref_en_(ref_en.end_),
      evt_en_(evt_en),
      total_len_(ref_en.length()) {}

ReadAln::ReadAln(const ReadAln &r)
    : ref_st_(r.ref_st_),
      evt_st_(r.evt_st_),
      ref_en_(r.ref_en_),
      evt_en_(r.evt_en_),
      total_len_(r.total_len_) {}


void ReadAln::update_next(ReadAln &new_loc) const {
    new_loc.ref_en_ = ref_en_;
    new_loc.evt_en_ = evt_en_;

    if (ref_st_.intersects(new_loc.ref_st_)) {

        if (new_loc.ref_st_.start_ < ref_st_.start_) {
            new_loc.total_len_ 
                = total_len_ + 
                  ref_st_.start_ - 
                  new_loc.ref_st_.start_;
        } else {
            new_loc.total_len_ = total_len_;
            new_loc.ref_st_.start_ = ref_st_.start_;
        }

    } else {
        new_loc.total_len_ += total_len_;
    }

}

Range ReadAln::ref_range() const {
    return Range(ref_st_.start_, ref_en_);
}

void ReadAln::print(std::ostream &out, bool newline = false, bool print_all = false) const {
    out << total_len_ << " ";

    out << ref_st_.start_;
    if (print_all ) {
        std::cout << ref_st_.end_ << ":";
    }

    out << "-" << ref_en_ << " " 
               << evt_st_ << "-" << evt_en_;

    if (newline)
        out << "\n";
}

bool ReadAln::is_valid() {
    return evt_st_ >= 0;
}

bool operator< (const ReadAln &r1, const ReadAln &r2) {
    if (r1.ref_st_.end_ != r2.ref_st_.end_)
        return r1.ref_st_.end_ < r2.ref_st_.end_;

    return r1.evt_st_ < r2.evt_st_;
}

SeedTracker::SeedTracker() {
    longest_seed = 0;
}

void SeedTracker::reset() {
    locations.clear();
    longest_seed = 0;
}

int SeedTracker::add_seeds(const std::vector<Result> &seeds) {
    int max_len = 0;
    for (int i = 0; i < seeds.size(); i++) {
        max_len = max(max_len, add_seed(seeds[i]));
    }
    return max_len;
}

int SeedTracker::add_seed(Result r) {
    ReadAln new_loc(r.ref_range_, r.read_range_.end_);

    auto loc = locations.lower_bound(new_loc),
         loc_match = locations.end();

    while (loc != locations.end()) {
        bool higher_sup = loc_match == locations.end() 
                       || loc_match->total_len_ < loc->total_len_,
             
             in_range = loc->ref_st_.end_ >= new_loc.ref_st_.end_
                        && loc->ref_st_.end_ - loc->evt_st_
                           <= new_loc.ref_st_.end_ - new_loc.evt_st_ + 6;

        if (higher_sup && in_range) {
            loc_match = loc;
        }

        loc++;
    }

    if (loc_match != locations.end()) {

        loc_match->update_next(new_loc);
        locations.erase(loc_match);
        locations.insert(new_loc);


        #ifdef DEBUG
        new_loc.print(true, true);
        #endif

    } else if (loc_match == locations.end()) {
        locations.insert(new_loc);

        #ifdef DEBUG
        new_loc.print(true, true);
        #endif
    }

    if (new_loc.total_len_ > longest_seed) {
        longest_seed = new_loc.total_len_;
    }

    //std::cerr << new_loc.ref_en_ << " " << r.ref_range_.length() << "\n";

    return new_loc.total_len_;
}

std::vector<ReadAln> SeedTracker::get_alignments(int min_len = 1) {
    std::map<int, ReadAln> sorted_alns;

    for (auto a = locations.begin(); a != locations.end(); a++) {
        if (a->total_len_ >= min_len) {
            sorted_alns[a->total_len_] = *a;
        }
    }

    std::vector<ReadAln> ret;
    for (auto a = sorted_alns.rbegin(); a != sorted_alns.rend(); a++) {
        ret.push_back(a->second);
    }

    return ret;
}

double SeedTracker::top_ratio(const std::vector<ReadAln> &alns) {
    const ReadAln &top = alns.front();

    Range top_ref = top.ref_range();

    int i = 1;
    for(; i < alns.size(); i++) {
        if (top_ref.get_recp_overlap(alns[i].ref_range()) < 0.9) {
            break;
        }
    }

    return double(top.total_len_) / double(alns[i].total_len_);
}

void SeedTracker::print(std::ostream &out, std::string strand, int max_out = 10) {

    std::vector<ReadAln> alns = get_alignments(1);

    //std::cout << strand << " " << top_ratio(alns) << "\n";

    if (alns.empty()) {
        return;
    }

    Range top_ref = alns[0].ref_range();
    double top_len = alns[0].total_len_;

    for (int i = 0; i < min(max_out, alns.size()); i++) {
        double overlap = top_ref.get_recp_overlap(alns[i].ref_range()),
               len_ratio = top_len / alns[i].total_len_;

        out << strand << " ";
        alns[i].print(out, false);
        out << "\t" << len_ratio << "\t" << overlap << "\n";
    }
}

