#ifndef INCL_NANOBWT
#define INCL_NANOBWT

#include "kmer_model.hpp"
#include <list>

class NanoFMI {
    public:

    class Range {
        
        public:
        int start_, end_;

        Range(NanoFMI &fmi);

        //Copy constructor
        Range(const Range &prev);

        Range(int start, int end);

        Range();

        Range& operator=(const Range &r);

        Range split_range(const Range &r);

        bool same_range(const Range &r) const;

        bool intersects(const Range &r) const;

        bool is_valid() const;

        friend bool operator< (const Range &q1, const Range &q2);
    };

    typedef struct Result {
        int qry_start, qry_end, ref_start, ref_end;
        float prob;
    } Result;

    NanoFMI(KmerModel &model, std::vector<mer_id> &mer_seq, int tally_dist);

    std::list<Range> get_neigbhors(Range range, std::list<mer_id> kmers) const;
    Range get_full_range(mer_id kmer) const;

    bool operator() (unsigned int rot1, unsigned int rot2);


    //private:
    //int tally_cp_dist(int i);
    //int get_tally(mer_id c, int i);
    std::list<int> get_tallies(std::list<mer_id> kmers, int loc) const;
    float get_stay_prob(Event e1, Event e2) const; //hmm

    KmerModel *model_;
    std::vector<mer_id> *mer_seq_;

    //Sparseness of suffix and tally arrays
    int tally_dist_;

    std::vector<mer_id> bwt_;                        //L array
    std::vector<int> mer_counts_, mer_f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<int> > mer_tally_;      //Rank tally
    mutable std::vector<int> mer_count_tmp_;


    void make_bwt();

    friend bool operator< (const NanoFMI::Range &q1, const NanoFMI::Range &q2);

    class SeedNode {
        public:
        int max_length_, stay_count_;
        mer_id kmer_;
        double prob_;


        std::list< std::pair<SeedNode *, bool> > parents_;
        std::list<SeedNode *> children_;
    
        //Source constructor
        SeedNode(mer_id kmer, double prob);

        //Child constructor
        SeedNode(SeedNode *parent, mer_id kmer, double prob, bool stay);

        //Creates invalid node
        SeedNode();

        //Copy constructor
        SeedNode(const SeedNode &s);

        ~SeedNode();

        //SeedNode& operator=(const SeedNode &s);

        bool is_valid();
        
        void invalidate(bool print=false);
        bool remove_child(SeedNode *child, bool invalidate, bool print=false);
        
        int add_child(SeedNode *child);


    };

    class SeedGraph {
        public:

        std::map<Range, SeedNode *> prev_nodes, next_nodes;
        std::list< std::list< SeedNode * > > sources_;
        std::list<double *> event_kmer_probs_;

        const NanoFMI &fmi_;
        const NormParams &norm_params_;
        unsigned int seed_length_, cur_event_;

        double event_prob_,
               seed_prob_,
               stay_prob_;
        
        SeedGraph(const NanoFMI &fmi, 
                  const NormParams &norm_params, 
                  int seed_len, int read_len,
                  double event_prob,
                  double seed_prob,
                  double stay_prob);

        std::vector<Result> add_event(Event e);
        void print_graph();

        SeedNode *add_child(Range &range, SeedNode &node); //copy update_ranges
        int add_sources(const Range &range, const SeedNode &node);


        std::vector<Result> pop_seeds(); //Big result gathering loop
                                    //Probably split into a few methods
        
        bool empty();
    };
    typedef std::pair<SeedNode *, bool> parent_ptr;
};

#endif

