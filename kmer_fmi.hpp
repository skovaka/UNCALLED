#ifndef INCL_NANOBWT
#define INCL_NANOBWT

#include "kmer_model.hpp"
#include <list>

class Range {
    
    public:
    int start_, end_;

    //Copy constructor
    Range(const Range &prev);

    Range(int start, int end);

    Range();

    Range& operator=(const Range &r);

    Range split_range(const Range &r);

    Range intersect(const Range &r) const;

    Range merge(const Range &r) const;

    double get_recp_overlap(const Range &r) const;

    bool same_range(const Range &r) const;

    bool intersects(const Range &r) const;

    bool is_valid() const;

    int length() const;

    friend bool operator< (const Range &q1, const Range &q2);
};

bool operator< (const Range &q1, const Range &q2);

class KmerFMI {
    public:


    KmerFMI(int alph_size, std::vector<Kmer> &mer_seq, int tally_dist);

    std::list<Range> get_neigbhors(Range range, std::list<Kmer> kmers) const;
    Range get_full_range(Kmer kmer) const;

    bool operator() (unsigned int rot1, unsigned int rot2);

    //private:
    std::list<int> get_tallies(std::list<Kmer> kmers, int loc) const;

    std::vector<Kmer> *mer_seq_;

    //Sparseness of suffix and tally arrays
    int tally_dist_;

    int alph_size_;
    std::vector<Kmer> bwt_;                        //L array
    std::vector<int> mer_counts_, mer_f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<int> > mer_tally_;      //Rank tally
    mutable std::vector<int> mer_count_tmp_;

    void make_bwt();

};

#endif

