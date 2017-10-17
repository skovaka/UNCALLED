#ifndef SEED_GRAPH_HPP
#define SEED_GRAPH_HPP

#include "nano_fmi.hpp"
#include "kmer_model.hpp"
#include <list>
#include <iostream>

//#define DEBUG_PROB

class Result {
    public:

    Result(int read_start, 
           int seed_len, 
           double prob, 
           int ref_start = 0, 
           int ref_end = 0);

    void set_ref_range(int start, int end);
    void print(std::ostream &out);

    Range read_range_, ref_range_;
    double seed_prob_;

    #ifdef DEBUG_PROB
    double min_evt_prob_;
    #endif
};

class AlnParams {
    public:
    const KmerModel &model_;

    int graph_elen_, anchor_rlen_, max_ignores_, max_skips_;

    double max_stay_frac_,
           min_anchor_evpr_,
           min_seed_pr_,
           min_stay_pr_,
           min_extend_evpr_;

    AlnParams(const KmerModel &model,
              int min_seed_nlen, 
              int anchor_nlen, 
              int max_ignores, 
              int max_skips,
              double max_stay_frac,
              double min_anchor_evpr,
              double min_extend_evpr,
              double min_seed_pr,
              double min_stay_pr);
    
    int nucl_to_events(int n);
    int get_graph_len(int seed_nlen);
};


class SeedGraph {

    class Node {
        public:

        enum Type { MATCH, STAY, SKIP, IGNORE };

        int length_, 
            stay_count_,
            skip_count_,
            ignore_count_;
        mer_id kmer_;
        double event_prob_, seed_prob_;

        #ifdef DEBUG_PROB
        double min_evt_prob_;
        double min_stay_prob_;
        #endif

        std::list< std::pair<Node *, Type> > parents_;
        std::list<Node *> children_;
    
        //Source constructor
        Node(mer_id kmer, double prob);

        //Child constructor
        Node(Node *parent, mer_id kmer, double prob, Type type);

        //Creates invalid node
        Node();

        //Copy constructor
        Node(const Node &s);

        void replace_info(const Node &node);
        void update_info();
        bool better_than(const Node *node); //TODO: this is a terrible name
        bool should_report(const AlnParams &params);

        int seed_len();
        int match_len();
        bool is_valid();

        
        void invalidate(std::vector<Node *> *old_nodes, bool delete_source);
        bool remove_child(Node *child);
        
        int add_child(Node *child);
    };
    typedef std::pair<Node *, Node::Type> parent_ptr;


    public:

    const NanoFMI &fmi_;

    AlnParams params_;
    NormParams norm_params_;
    std::string label_;

    std::map<Range, Node *> next_nodes_;
    std::map<Node *, Range> prev_nodes_;
    std::list< std::list< Node * > > sources_;
    std::list<double *> event_kmer_probs_;
    std::vector<Node *> old_nodes_;
    
    unsigned int cur_event_;
    Event prev_event_;

    SeedGraph(const NanoFMI &fmi, 
              const AlnParams &aln_params,
              const std::string &label);

    ~SeedGraph();

    void new_read(int read_len, const NormParams &params);
    void reset();

    std::vector<Result> add_event(Event e, std::ostream &out);
    void print_graph(bool verbose);

    Node *add_child(Range &range, Node &node); //copy update_ranges
    int add_sources(const Range &range, const Node &node);

    std::vector<Result> pop_seeds(std::ostream &out); //Big result gathering loop
                                     //Probably split into a few methods
    bool empty();
};


#endif

