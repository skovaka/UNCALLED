#include "seed_tracker.hpp"
#include <iostream>
#include <set>

int main(int argc, char **argv) {
    
    std::string prev_strand, strand, evt_str, ref_str;
    double prob;

    int i, evt_st, evt_en, ref_st, ref_en;

    #ifdef DEBUG_PROB
    double min_evt_prob_;
    #endif

    SeedTracker fwd_tracker, rev_tracker;

    prev_strand = "";

    while (std::cin) {

        std::cin >> strand;
        
        if (strand[0] == '=') {
            getline(std::cin, strand);
            strand = prev_strand;
            continue;
        }

        std::cin >> evt_str >> ref_str >> prob;

        #ifdef DEBUG_PROB
        std::cin >> min_evt_prob_;
        #endif

        i = evt_str.find('-');
        evt_st = atoi(evt_str.substr(0, i).c_str());
        evt_en = atoi(evt_str.substr(i + 1).c_str());

        i = ref_str.find('-');
        ref_st = atoi(ref_str.substr(0, i).c_str());
        ref_en = atoi(ref_str.substr(i + 1).c_str());

        Result r(evt_st, evt_en-evt_en, prob, ref_st, ref_en);

        if (strand == "fwd") {
            fwd_tracker.add_seed(r);
        } else {
            rev_tracker.add_seed(r);
        }

        #ifdef DEBUG_PROB
        std::cerr << ref_en << "\t" << min_evt_prob_ << "\n";
        #endif

        if (strand[0] != '=') {
            prev_strand = strand;
        }

    }

    std::string s1 = "fwd", s2 = "rev";
    fwd_tracker.print(s1);
    rev_tracker.print(s2);
}
