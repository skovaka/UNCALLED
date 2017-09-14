#ifndef SEED_TRACKER_HPP
#define SEED_TRACKER_HPP

#include "seed_graph.hpp"
#include <set>

class SeedTracker {
    
    class RefLoc {
        public:
        int ref_st_, evt_st_, sup_,
            ref_en_, evt_en_;

        RefLoc(int ref_en, int evt_en);
        RefLoc(const RefLoc &r);

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
