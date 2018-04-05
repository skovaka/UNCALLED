#ifndef SEED_GRAPH_HPP
#define SEED_GRAPH_HPP

#include "fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"
#include <list>
#include <iostream>

//#define DEBUG_NODES
//#define DEBUG_PROB

class Result {
    public:

    Result(unsigned int read_start, 
           unsigned int seed_len, 
           double prob, 
           unsigned int ref_start = 0, 
           unsigned int ref_end = 0);

    void set_ref_range(unsigned int start, unsigned int end);
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

    size_t graph_elen_, anchor_rlen_, max_ignores_, max_skips_, max_consec_stay_;

    double max_stay_frac_,
           min_anchor_evpr_,
           min_seed_pr_,
           min_stay_pr_;
           //min_extend_evpr_;

    std::vector<unsigned int> expr_lengths_;
    std::vector<double> expr_probs_;

    AlnParams(const KmerModel &model,
              unsigned int min_seed_nlen, 
              unsigned int anchor_nlen, 
              unsigned int max_ignores, 
              unsigned int max_skips,
              unsigned int max_consec_stay,
              double max_stay_frac,
              double min_anchor_evpr,
              //double min_extend_evpr,
              std::vector<unsigned int> expr_lengths,
              std::vector<double> expr_probs,
              double min_seed_pr,
              double min_stay_pr);
    
    unsigned int nucl_to_events(unsigned int n);
    unsigned int get_graph_len(unsigned int seed_nlen);
};


class AlignmentForest {

    class Node {
        public:

        enum Type { MATCH, STAY, SKIP, IGNORE };

        size_t length_, 
               stay_count_,
               skip_count_,
               ignore_count_,
               consec_stays_;

        unsigned char child_count_;

        mer_id kmer_;

        double event_prob_, seed_prob_;
        Node::Type type_;

        Node *parent_;
        Node **children_;
    
        //Source constructor
        Node(mer_id kmer, double prob);

        //Child constructor
        Node(Node *parent, mer_id kmer, double prob, Type type);

        //Creates invalid node
        Node();

        //Copy constructor
        Node(const Node &s);
        ~Node();

        void update_info();
        bool better_than(const Node *node); //TODO: this is a terrible name
        bool should_report(const AlnParams &params);

        size_t seed_len();
        size_t match_len();
        bool is_valid();
        double mean_prob() const;
        double next_mean_prob(double next_prob) const;

        void invalidate(std::vector<Node *> *old_nodes, bool delete_source);
        bool remove_child(Node *child);
        size_t add_child(Node *child);
        void print() const;
    };

    //class SingleNode : Node {
    //}
    //typedef std::pair<Node *, Node::Type> parent_ptr;


    public:

    const FMI &fmi_;
    Range *kmer_ranges_;
    std::map<Range, size_t> range_kmers_;

    AlnParams params_;
    std::string label_;


    std::map<Range, Node *> next_nodes_;
    std::map<Node *, Range> prev_nodes_;
    std::list< std::list< Node * > > sources_;
    std::list<double *> event_kmer_probs_;
    std::vector<Node *> old_nodes_;

    Timer timer;
    
    unsigned int cur_event_;
    Event prev_event_;

    AlignmentForest(const FMI &fmi, 
              const AlnParams &aln_params,
              const std::string &label);

    ~AlignmentForest();

    void new_read(size_t read_len);
    void reset();

    std::vector<Result> add_event(Event e, std::ostream &out);
    void print_graph(bool verbose);

    Node *add_child(Range &range, Node &node); //copy update_ranges
    size_t add_sources(const Range &range, const Node &node);

    std::vector<Result> pop_seeds(std::ostream &out); //Big result gathering loop
                                     //Probably split into a few methods
    bool empty();
};


#endif

