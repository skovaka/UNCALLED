#ifndef SEED_GRAPH_HPP
#define SEED_GRAPH_HPP

#include "nano_fmi.hpp"
#include "kmer_model.hpp"
#include <list>

typedef struct Result {
    int qry_start, qry_end, ref_start, ref_end;
    float prob;
} Result;



class SeedGraph {
    class Node {
        public:
        int max_length_, stay_count_;
        mer_id kmer_;
        double prob_;


        std::list< std::pair<Node *, bool> > parents_;
        std::list<Node *> children_;
    
        //Source constructor
        Node(mer_id kmer, double prob);

        //Child constructor
        Node(Node *parent, mer_id kmer, double prob, bool stay);

        //Creates invalid node
        Node();

        //Copy constructor
        Node(const Node &s);

        ~Node();

        //Node& operator=(const Node &s);

        bool is_valid();
        
        void invalidate(bool print=false);
        bool remove_child(Node *child, bool invalidate, bool print=false);
        
        int add_child(Node *child);
    };

    public:

    std::map<Range, Node *> prev_nodes, next_nodes;
    std::list< std::list< Node * > > sources_;
    std::list<double *> event_kmer_probs_;
    
    const KmerModel &model_;
    const NanoFMI &fmi_;
    const NormParams &norm_params_;
    unsigned int seed_length_, cur_event_, max_stays_;

    double event_prob_,
           seed_prob_,
           stay_prob_;
    
    SeedGraph(const KmerModel &model,
              const NanoFMI &fmi, 
              const NormParams &norm_params, 
              int seed_len, int read_len,
              double event_prob,
              double seed_prob,
              double stay_prob,
              double stay_frac);

    std::vector<Result> add_event(Event e);
    void print_graph();

    Node *add_child(Range &range, Node &node); //copy update_ranges
    int add_sources(const Range &range, const Node &node);


    std::vector<Result> pop_seeds(); //Big result gathering loop
                                //Probably split into a few methods
    
    bool empty();



    typedef std::pair<Node *, bool> parent_ptr;
};


#endif

