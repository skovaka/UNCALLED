#ifndef INCL_BASEFMI
#define INCL_BASEFMI

#include <list>
#include <vector>
#include <string>
#include "fmi.hpp"
#include "basepairs.hpp"
#include "range.hpp"


class BaseFMI : public FMI {
    public:

    BaseFMI(std::string seq, unsigned int tally_dist);
    BaseFMI();
    BaseFMI(std::ifstream &infile, unsigned int tally_dist);

    void construct(const std::string &seq);

    void save(const std::string &filename); 

    Range get_neighbor(Range range, base_t base) const;

    Range get_full_range(base_t base) const;

    Range get_kmer_range(const std::string &seq) const;

    size_t sa(size_t i) const;

    size_t size() const;

    bool operator() (unsigned int rot1, unsigned int rot2);

    //private:
    unsigned int get_tally(base_t bases, unsigned int loc) const;

    std::string *seq_;

    //Sparseness of suffix and tally arrays
    unsigned int tally_dist_;

    std::vector<base_t> bwt_;                        //L array
    std::vector<unsigned int> counts_, f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<unsigned int> > tally_;      //Rank tally

    void make_bwt();

};

#endif
