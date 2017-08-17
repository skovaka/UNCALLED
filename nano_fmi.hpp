#ifndef INCL_NANOBWT
#define INCL_NANOBWT

#include "kmer_model.hpp"
#include <list>

class NanoFMI 
{
    public:
    typedef struct Result {
        int qry_start, qry_end, ref_start, ref_end;
        float prob;
    } Result;

    NanoFMI(KmerModel &model, std::vector<mer_id> &mer_seq, int tally_dist);

    bool operator() (unsigned int rot1, unsigned int rot2);

    std::vector<Result> lf_map(std::vector<Event> &events, int seed_end, 
               int match_len, NormParams norm_params);




    private:
    int tally_cp_dist(int i);
    int get_tally(mer_id c, int i);
    float get_stay_prob(Event e1, Event e2); //hmm

    KmerModel *model_;
    std::vector<mer_id> *mer_seq_;



    //Sparseness of suffix and tally arrays
    int tally_dist_;

    std::vector<mer_id> bwt_;                        //L array
    std::vector<int> mer_counts_, mer_f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<int> > mer_tally_;      //Rank tally


    void make_bwt();

    //typedef struct RangeMetadata {
    //     int match_len_;              
    //     double prob_;
    //     bool traversed_;
    //     Range const * parent_, self_;
    //     std::list<Range const *> children_;
    //} RangeMetadata;
    
    class Range {
        //private:


        public:
        NanoFMI &fmi_;
        mer_id k_id_;
        int event_, start_, end_;
        mutable int match_len_, total_len_;
        mutable double prob_;
        mutable bool stay_;
        mutable Range const * parent_;
        mutable int child_count_;

        Range(NanoFMI &fmi);

        //Copy constructor
        Range(const Range &prev);

        //Initial match constructor
        Range(NanoFMI &fmi, mer_id k_id, int event, double prob);

        //"next" constructor
        Range(const Range &prev, mer_id k_id, double prob);

        //"stay" constructor
        Range(const Range &prev, double prob);

        ~Range();

        Range& operator=(const Range&);

        unsigned int id() const;

        bool add_results(std::vector<Result> &results, 
                         int range_end, double min_prob) const;

        bool same_range(const Range &q) const;

        bool is_valid() const;

        bool is_match() const;

        void print_info() const;

        bool intersects(const Range &q) const;

        Range split_range(const Range &q);

        friend bool operator< (const Range &q1, const Range &q2);


    };
    friend bool operator< (const NanoFMI::Range &q1, const NanoFMI::Range &q2);
    friend bool check_ranges(const NanoFMI::Range &q1, const NanoFMI::Range &q2);

    int update_ranges(std::set<Range> &next_ranges, 
                      const Range &nr, bool split_all);

    typedef std::vector< std::set<Range> > EventQueries;
};

#endif

