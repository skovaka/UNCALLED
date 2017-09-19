#include "seed_tracker.hpp"
#include "nano_fmi.hpp"
#include <iostream>
#include <set>

SeedTracker::RefLoc::RefLoc(int ref_en, int evt_en, int len)
    : ref_st_(ref_en),
      evt_st_(evt_en),
      ref_en_(ref_en),
      evt_en_(evt_en),
      prev_len_(len), 
      total_len_(len) {}

SeedTracker::RefLoc::RefLoc(const RefLoc &r)
    : ref_st_(r.ref_st_),
      evt_st_(r.evt_st_),
      ref_en_(r.ref_en_),
      evt_en_(r.evt_en_),
      prev_len_(r.prev_len_), 
      total_len_(r.total_len_) {}

int SeedTracker::RefLoc::ref_start_base() const {
    return ref_st_ - prev_len_;
}

void SeedTracker::RefLoc::update_next(RefLoc &new_loc) const {
    new_loc.ref_en_ = ref_en_;
    new_loc.evt_en_ = evt_en_;

    std::cerr << new_loc.prev_len_ << "\n";

    if (new_loc.ref_start_base() < ref_start_base()) {
        if (new_loc.ref_st_ >= ref_start_base()) {

            new_loc.total_len_ 
                = total_len_ + 
                  ref_start_base() - 
                  new_loc.ref_start_base();

        } else {
            new_loc.total_len_ += total_len_;
        }
    } else {
        new_loc.total_len_ = total_len_;
    }

}

void SeedTracker::RefLoc::print() const {
    std::cout << total_len_ << " " 
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
    RefLoc new_loc(r.ref_range_.end_, r.read_range_.end_, r.ref_range_.length());

    auto loc = locations.lower_bound(new_loc),
         loc_match = locations.end();

    while (loc != locations.end()) {
        bool higher_sup = loc_match == locations.end() 
                       || loc_match->total_len_ < loc->total_len_,
             
             in_range = loc->ref_st_ >= new_loc.ref_st_
                        && loc->ref_st_ - loc->evt_st_
                           <= new_loc.ref_st_ - new_loc.evt_st_ + 2;

        if (higher_sup && in_range) {
            loc_match = loc;
        }

        loc++;
    }


    if (loc_match != locations.end() && 
        new_loc.evt_st_ != loc_match->evt_st_) {

        loc_match->update_next(new_loc);
        locations.erase(loc_match);
        locations.insert(new_loc);
        //new_loc.print();

    } else if (loc_match == locations.end()) {
        locations.insert(new_loc);
        //new_loc.print();
    }

    //std::cerr << new_loc.ref_en_ << " " << r.ref_range_.length() << "\n";

    return new_loc.ref_en_;
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

    #ifdef DEBUG_PROB
    double min_evt_prob_;
    #endif

    SeedTracker tracker;

    prev_strand = "";

    while (std::cin) {
        std::cin >> strand;
        
        if (strand[0] == '=') {
            getline(std::cin, strand);
            continue;
        }

        std::cin >> evt_str >> ref_str >> prob;

        #ifdef DEBUG_PROB
        std::cin >> min_evt_prob_;
        #endif

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

        int ref_en = tracker.add_seed(r);

        #ifdef DEBUG_PROB
        std::cerr << ref_en << "\t" << min_evt_prob_ << "\n";
        #endif

        prev_strand = strand;
    }

    tracker.print(prev_strand);

}
