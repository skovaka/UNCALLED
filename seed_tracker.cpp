#include "seed_tracker.hpp"
#include "nano_fmi.hpp"
#include <iostream>
#include <set>


SeedTracker::RefLoc::RefLoc(int ref_en, int evt_en)
    : ref_st_(ref_en),
      evt_st_(evt_en),
      sup_(1),
      ref_en_(ref_en),
      evt_en_(evt_en) {}

SeedTracker::RefLoc::RefLoc(const RefLoc &r)
    : ref_st_(r.ref_st_),
      evt_st_(r.evt_st_),
      sup_(r.sup_),
      ref_en_(r.ref_en_),
      evt_en_(r.evt_en_) {}

void SeedTracker::RefLoc::update_next(RefLoc &new_loc) const {
    new_loc.ref_en_ = ref_en_;
    new_loc.evt_en_ = evt_en_;
    new_loc.sup_ = sup_ + 1;
}

void SeedTracker::RefLoc::print() const {
    std::cout << sup_ << " " 
              << ref_st_ << "-" << ref_en_ << " " 
              << evt_st_ << "-" << evt_en_ << "\n";
}

bool operator< (const SeedTracker::RefLoc &r1, const SeedTracker::RefLoc &r2) {
    if (r1.ref_st_ != r2.ref_st_)
        return r1.ref_st_ < r2.ref_st_;

    return r1.evt_st_ < r2.evt_st_;
}

SeedTracker::SeedTracker() {

}

int SeedTracker::add_seed(Result r) {
    RefLoc new_loc(r.ref_range_.end_, r.read_range_.end_);

    auto loc = locations.lower_bound(new_loc),
         loc_match = locations.end();


    while (loc != locations.end()) {
        bool higher_sup = loc_match == locations.end() 
                          || loc_match->sup_ < loc->sup_,
             
             in_range = loc->ref_st_ >= new_loc.ref_st_
                        && loc->ref_st_ - loc->evt_st_ 
                           <= new_loc.ref_st_ - new_loc.evt_st_;

        if (higher_sup && in_range) {
            std::cout << (loc->ref_st_ < new_loc.ref_st_) << "\n";
            loc_match = loc;
        }

        loc++;
    }


    if (loc_match != locations.end() && 
        new_loc.evt_st_ != loc_match->evt_st_) {

        loc_match->update_next(new_loc);
        locations.erase(loc_match);
        locations.insert(new_loc);

    } else if (loc_match == locations.end()) {
        locations.insert(new_loc);
    }

    return new_loc.sup_;
}

void SeedTracker::print(std::string &strand) {
    for (auto l = locations.begin(); l != locations.end(); l++) {
        std::cout << strand << " ";
        l->print();
    }
}

int main(int argc, char **argv) {
    
    std::string prev_strand, strand, evt_str, ref_str;
    double prob;

    int i, evt_st, evt_en, ref_st, ref_en;

    SeedTracker tracker;

    prev_strand = "";

    while (std::cin) {
        std::cin >> strand;
        
        if (strand[0] == '=') {
            getline(std::cin, strand);
            continue;
        }

        std::cin >> evt_str >> ref_str >> prob;

        if (!prev_strand.empty() && prev_strand != strand) {
            tracker.print(prev_strand);     
            tracker = SeedTracker();
        }

        i = evt_str.find('-');
        evt_st = atoi(evt_str.substr(0, i).c_str());
        evt_en = atoi(evt_str.substr(i + 1).c_str());

        i = ref_str.find('-');
        ref_st = atoi(ref_str.substr(0, i).c_str());
        ref_en = atoi(ref_str.substr(i + 1).c_str());

        Result r(evt_st, evt_en-evt_en, prob, ref_st, ref_en);

        tracker.add_seed(r);
        prev_strand = strand;
    }

    tracker.print(prev_strand);

}
