#ifndef SEED_TRACKER_HPP
#define SEED_TRACKER_HPP

#include "seed_graph.hpp"
#include <set>

class SeedTracker {
    
    class RefLoc {
        public:
        int ref_st_, evt_st_,
            ref_en_, evt_en_,
            prev_len_, total_len_;

        RefLoc(int ref_en, int evt_en, int len);
        RefLoc(const RefLoc &r);
        int ref_start_base() const;
        void update_next(RefLoc &new_loc) const;
        void print() const;

        friend bool operator< (const RefLoc &q1, const RefLoc &q2);
    };

    public:

    friend bool operator< (const RefLoc &q1, const RefLoc &q2);

    std::set<RefLoc> locations;

    SeedTracker();
    int add_seed(Result r);
    void print(std::string &strand);
};

#endif
