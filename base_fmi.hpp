#ifndef INCL_BASEFMI
#define INCL_BASEFMI

#include <list>
#include <vector>
#include <string>

#define ALPH_SIZE 4

class Range {
    
    public:
    unsigned int start_, end_;

    //Copy constructor
    Range(const Range &prev);

    Range(unsigned int start, unsigned int end);

    Range();

    Range& operator=(const Range &r);

    Range split_range(const Range &r);

    Range intersect(const Range &r) const;

    Range merge(const Range &r) const;

    double get_recp_overlap(const Range &r) const;

    bool same_range(const Range &r) const;

    bool intersects(const Range &r) const;

    bool is_valid() const;

    unsigned int length() const;

    friend bool operator< (const Range &q1, const Range &q2);
};

bool operator< (const Range &q1, const Range &q2);

class BaseFMI {
    public:

    BaseFMI(std::string seq, unsigned int tally_dist);
    BaseFMI();
    BaseFMI(std::ifstream &infile, unsigned int tally_dist);
    void save(std::string filename); 

    Range get_neighbor(Range range, char base) const;

    Range get_full_range(char base) const;

    Range get_kmer_range(const std::string &seq) const;

    bool operator() (unsigned int rot1, unsigned int rot2);

    //private:
    unsigned int get_tally(char bases, unsigned int loc) const;

    std::string *seq_;

    bool loaded_;

    //Sparseness of suffix and tally arrays
    unsigned int tally_dist_;

    std::string bwt_;                        //L array
    std::vector<unsigned int> counts_, f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<unsigned int> > tally_;      //Rank tally

    void make_bwt();

};

#endif
