#ifndef SEED_TRACKER_HPP
#define SEED_TRACKER_HPP

#include "seed_graph.hpp"
#include <set>
#include <vector>

unsigned int max(unsigned int a, unsigned int b);
unsigned int min(unsigned int a, unsigned int b);

class ReadAln {
    public:
    Range ref_st_;
    unsigned int evt_st_,
        ref_en_, evt_en_,
        total_len_;

    ReadAln(Range ref_en, unsigned int evt_en);
    ReadAln(const ReadAln &r);
    ReadAln();
    unsigned int ref_start_base() const;
    void update_next(ReadAln &new_loc) const;
    void print(std::ostream &out, bool newline, bool print_all) const;
    Range ref_range() const;
    bool is_valid();

    friend bool operator< (const ReadAln &q1, const ReadAln &q2);
};

bool operator< (const ReadAln &q1, const ReadAln &q2);

class SeedTracker {

    public:

    std::set<ReadAln> locations;
    ReadAln *longest_seed;
    unsigned int min_revents_;

    SeedTracker();
    SeedTracker(unsigned int min_revents);

    ReadAln add_seed(Result seed);
    ReadAln add_seeds(const std::vector<Result> &seeds);

    void reset();

    std::vector<ReadAln> get_alignments(unsigned int min_len);

    //static double top_ratio(int min_len);
    bool check_ratio(const ReadAln &aln, double ratio);

    void print(std::ostream &out, std::string strand, size_t max_out);
};


#endif
