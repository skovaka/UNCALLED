#ifndef INCL_BASEFMI
#define INCL_BASEFMI

#include <list>
#include <vector>
#include <string>

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

class BaseFMI {
    public:

    BaseFMI(std::string seq, int tally_dist);

    std::list<Range> get_neigbhors(Range range, std::list<char> bases) const;
    Range get_neighbor(Range range, char base) const;

    Range get_full_range(char base) const;

    Range get_kmer_range(const std::string &seq) const;

    bool operator() (unsigned int rot1, unsigned int rot2);

    //private:
    std::list<int> get_tallies(std::list<char> bases, int loc) const;
    int get_tally(char bases, int loc) const;

    std::string *seq_;

    //Sparseness of suffix and tally arrays
    int tally_dist_;

    int alph_size_;
    std::string bwt_;                        //L array
    std::vector<int> counts_, f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<int> > tally_;      //Rank tally
    mutable std::vector<int> count_tmp_;

    void make_bwt();

};

#endif
