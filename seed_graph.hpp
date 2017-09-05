#ifndef SEED_GRAPH_HPP
#define SEED_GRAPH_HPP

#include "nano_fmi.hpp"
#include "kmer_model.hpp"
#include <list>

class Result {
    public:

    Result(int read_start, int seed_len, double prob, int ref_start = 0, int ref_end = 0);

    void set_ref_range(int start, int end);
    void print();

    Range read_range_, ref_range_;
    double prob_;
};



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

        bool is_valid();
        
        void invalidate(std::vector<Node *> *old_nodes);
        bool remove_child(Node *child, std::vector<Node *> *old_nodes);
        
        int add_child(Node *child);
    };

    public:

    std::map<Range, Node *> next_nodes_;
    std::map<Node *, Range> prev_nodes_;
    std::list< std::list< Node * > > sources_;
    std::list<double *> event_kmer_probs_;
    std::vector<Node *> old_nodes_;
    
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

    ~SeedGraph();

    std::vector<Result> add_event(Event e);
    void print_graph(bool verbose);


    Node *add_child(Range &range, Node &node); //copy update_ranges
    int add_sources(const Range &range, const Node &node);


    std::vector<Result> pop_seeds(); //Big result gathering loop
                                //Probably split into a few methods
    
    bool empty();



    typedef std::pair<Node *, bool> parent_ptr;
};


#endif

