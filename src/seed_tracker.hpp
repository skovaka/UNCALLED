#ifndef SEED_TRACKER_HPP
#define SEED_TRACKER_HPP

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include "util.hpp"
#include "range.hpp"

class SeedGroup {
    public:

    u64 ref_st_;
    Range ref_en_;
    u32 evt_st_,
        evt_en_,
        total_len_;


    SeedGroup(Range ref_st, u32 evt_st);
    SeedGroup(const SeedGroup &r);
    SeedGroup();
    u64 ref_start_base() const;
    u8 update(SeedGroup &new_aln);
    void print(std::ostream &out, bool newline, bool print_all) const;
    Range ref_range() const;
    bool is_valid();

    friend bool operator< (const SeedGroup &q1, const SeedGroup &q2);
    friend std::ostream &operator<< (std::ostream &out, const SeedGroup &a);
};

const SeedGroup NULL_ALN = SeedGroup();

//static const SeedGroup NULL_ALN;

bool operator< (const SeedGroup &q1, const SeedGroup &q2);
std::ostream &operator<< (std::ostream &out, const SeedGroup &a);

class SeedTracker {

    public:

    std::set<SeedGroup> alignments_;
    std::multiset<u32> all_lens_;
    const u64 ref_len_;
    const float mean_thresh_, top_thresh_;
    const u8 min_aln_len_;

    float max_len_, len_sum_;

    SeedTracker(u64 ref_len, float mean_thresh, float top_thresh, u8 min_aln_len, u8 win_len);

    SeedGroup add_seed(SeedGroup sg);
    SeedGroup add_seeds(const std::vector<SeedGroup> &seeds);

    void reset();

    std::vector<SeedGroup> get_alignments(u8 min_len);

    //static double top_ratio(int min_len);
    bool check_ratio(const SeedGroup &aln, double ratio);

    void print(std::ostream &out, u16 max_out);
};


#endif
