#ifndef SEED_TRACKER_HPP
#define SEED_TRACKER_HPP

#include "seed_graph.hpp"
#include <set>
#include <vector>

int max(int a, int b);
int min(int a, int b);

class ReadAln {
    public:
    Range ref_st_;
    int evt_st_,
        ref_en_, evt_en_,
        total_len_;

    ReadAln(Range ref_en, int evt_en);
    ReadAln(const ReadAln &r);
    ReadAln();
    int ref_start_base() const;
    void update_next(ReadAln &new_loc) const;
    void print(std::ostream &out, bool newline, bool print_all) const;
    Range ref_range() const;

    friend bool operator< (const ReadAln &q1, const ReadAln &q2);
};

bool operator< (const ReadAln &q1, const ReadAln &q2);

class SeedTracker {

    public:


    std::set<ReadAln> locations;
    int longest_seed;

    SeedTracker();
    int add_seed(Result seed);
    int add_seeds(const std::vector<Result> &seeds);

    void reset();

    std::vector<ReadAln> get_alignments(int min_len);

    static double top_ratio(const std::vector<ReadAln> &alns);

    void print(std::ostream &out, std::string strand, int max_out);
};


#endif
