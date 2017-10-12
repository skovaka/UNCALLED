#include "seed_tracker.hpp"
#include "nano_fmi.hpp"
#include <iostream>
#include <set>

//#define DEBUG

ReadAln::ReadAln() {
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

void ReadAln::print(bool print_all = false) const {
    std::cout << total_len_ << " ";

    if (print_all ) {
        std::cout << ref_st_.start_ << ":";
    }

    std::cout << ref_st_.end_ << "-" << ref_en_ << " " 
              << evt_st_ << "-" << evt_en_ << "\n";
}

bool operator< (const ReadAln &r1, const ReadAln &r2) {
    if (r1.ref_st_.end_ != r2.ref_st_.end_)
        return r1.ref_st_.end_ < r2.ref_st_.end_;

    return r1.evt_st_ < r2.evt_st_;
}

SeedTracker::SeedTracker() {
}

int SeedTracker::add_seeds(const std::vector<Result> &seeds) {
    for (int i = 0; i < seeds.size(); i++) {
        add_seed(seeds[i]);
    }
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

        //std::cout << (loc->ref_st_.end_ >= new_loc.ref_st_.end_) << " " 
        //          << loc->ref_st_.end_ - loc->evt_st_ << " "
        //          << new_loc.ref_st_.end_ - new_loc.evt_st_ << " "
        //          << higher_sup << "\n";
        //loc->print();

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

        #ifdef DEBUG
        new_loc.print(true);
        #endif

    } else if (loc_match == locations.end()) {
        locations.insert(new_loc);

        #ifdef DEBUG
        new_loc.print(true);
        #endif
    }

    //std::cerr << new_loc.ref_en_ << " " << r.ref_range_.length() << "\n";

    

    return new_loc.ref_en_;
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

void SeedTracker::print(std::string &strand) {
    std::vector<ReadAln> alns = get_alignments(1);
    for (int i = 0; i < alns.size(); i++) {
        std::cout << strand << "\t";
        alns[i].print();
    }
}

//int main(int argc, char **argv) {
//    
//    std::string prev_strand, strand, evt_str, ref_str;
//    double prob;
//
//    int i, evt_st, evt_en, ref_st, ref_en;
//
//    #ifdef DEBUG_PROB
//    double min_evt_prob_;
//    #endif
//
//    SeedTracker tracker;
//
//    prev_strand = "";
//
//    while (std::cin) {
//
//        std::cin >> strand;
//        
//        if (strand[0] == '=') {
//            getline(std::cin, strand);
//            strand = prev_strand;
//            continue;
//        }
//
//        std::cin >> evt_str >> ref_str >> prob;
//
//        #ifdef DEBUG_PROB
//        std::cin >> min_evt_prob_;
//        #endif
//
//        if (!prev_strand.empty() && prev_strand != strand) {
//            tracker.print(prev_strand);     
//            tracker = SeedTracker();
//        }
//
//        i = evt_str.find('-');
//        evt_st = atoi(evt_str.substr(0, i).c_str());
//        evt_en = atoi(evt_str.substr(i + 1).c_str());
//
//        i = ref_str.find('-');
//        ref_st = atoi(ref_str.substr(0, i).c_str());
//        ref_en = atoi(ref_str.substr(i + 1).c_str());
//
//        Result r(evt_st, evt_en-evt_en, prob, ref_st, ref_en);
//
//        int ref_en = tracker.add_seed(r);
//
//        #ifdef DEBUG_PROB
//        std::cerr << ref_en << "\t" << min_evt_prob_ << "\n";
//        #endif
//
//        if (strand[0] != '=') {
//            prev_strand = strand;
//        }
//
//    }
//
//    tracker.print(prev_strand);
//
//}
