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

    //typedef struct QueryMetadata {
    //     int match_len_;              
    //     double prob_;
    //     bool traversed_;
    //     Query const * parent_, self_;
    //     std::list<Query const *> children_;
    //} QueryMetadata;
    
    class Query {
        //private:


        public:
        NanoFMI &fmi_;
        mer_id k_id_;
        int start_, end_;
        mutable int match_len_;
        mutable double prob_;
        mutable bool traversed_;
        mutable Query const * parent_;
        mutable std::list<Query const *> children_;

        Query(NanoFMI &fmi);

        //Copy constructor
        Query(const Query &prev);

        //Initial match constructor
        Query(NanoFMI &fmi, mer_id k_id, double prob);

        //"next" constructor
        Query(const Query &prev, mer_id k_id, double prob);

        //"stay" constructor
        Query(const Query &prev, double prob);

        ~Query();

        Query& operator=(const Query&);

        unsigned int id() const;

        bool add_results(std::vector<Result> &results, 
                         int query_end, double min_prob) const;

        bool same_range(const Query &q) const;

        bool is_valid() const;

        bool is_match() const;

        void print_info() const;

        bool intersects(const Query &q) const;

        Query split_query(Query &q);

        friend bool operator< (const Query &q1, const Query &q2);


    };
    friend bool operator< (const NanoFMI::Query &q1, const NanoFMI::Query &q2);
    friend bool check_queries(const NanoFMI::Query &q1, const NanoFMI::Query &q2);

    int update_queries(std::set<Query> &queries, Query &q, bool first);

    typedef std::vector< std::set<Query> > EventQueries;
};

#endif

