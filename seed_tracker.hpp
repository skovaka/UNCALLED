#ifndef SEED_TRACKER_HPP
#define SEED_TRACKER_HPP

#include "aligner.hpp"
#include <set>
#include <vector>

unsigned int max(unsigned int a, unsigned int b);
unsigned int min(unsigned int a, unsigned int b);

class ReadAln {
    public:
    Range ref_en_;
    unsigned int evt_st_,
                 ref_st_, 
                 evt_en_,
                 total_len_;

    float prob_sum_;

    u8 segments_;

    ReadAln(Range ref_st, unsigned int evt_st, float prob);
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
    unsigned int read_length_, min_revents_;

    SeedTracker(unsigned int read_length);

    ReadAln add_seed(Result seed);
    ReadAln add_seeds(const std::vector<Result> &seeds);

    void reset();

    std::vector<ReadAln> get_alignments(unsigned int min_len);

    //static double top_ratio(int min_len);
    bool check_ratio(const ReadAln &aln, double ratio);

    void print(std::ostream &out, std::string strand, size_t max_out);
};


#endif
