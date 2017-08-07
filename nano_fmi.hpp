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
    
    class Query {
        public:

        //Copy constructor
        Query(const Query &prev);

        //Initial match constructor
        Query(NanoFMI &fmi, mer_id k_id, float prob);

        //"next" constructor
        Query(const Query &prev, mer_id k_id, double prob);

        //"stay" constructor
        Query(const Query &prev, float prob);

        bool add_results(std::vector<Result> &results, 
                         int query_end, double min_prob) const;

        bool same_range(const Query &q) const;

        bool is_valid() const;

        int match_len() const;

        float avg_prob() const;

        friend bool operator< (const Query &q1, const Query &q2);

        private:

        NanoFMI &fmi_;
        int start_, end_, match_len_, stays_;
        float prob_sum_;

    };
    friend bool operator< (const NanoFMI::Query &q1, const NanoFMI::Query &q2);

    typedef std::vector< std::set<Query> > EventQueries;
};

#endif

