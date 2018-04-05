#ifndef GRAPH_ALIGNER_HPP
#define GRAPH_ALIGNER_HPP

#include "aligner.hpp"
#include "sdsl_fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"
#include <list>
#include <iostream>

//#define DEBUG_NODES
//#define DEBUG_PROB

class GraphAligner : Aligner {

    class Node {
        public:

        enum Type { MATCH, STAY, SKIP, IGNORE };

        size_t length_, 
               stay_count_,
               skip_count_,
               ignore_count_,
               consec_stays_;

        Kmer kmer_;
        double event_prob_, seed_prob_;

        #ifdef DEBUG_NODES
        size_t full_length_;
        #endif

        #ifdef DEBUG_PROB
        double min_evt_prob_;
        double min_stay_prob_;
        #endif

        std::list< std::pair<Node *, Type> > parents_;
        std::list<Node *> children_;
    
        //Source constructor
        Node(Kmer kmer, double prob);

        //Child constructor
        Node(Node *parent, Kmer kmer, double prob, Type type);

        //Creates invalid node
        Node();

        //Copy constructor
        Node(const Node &s);


        void replace_info(const Node &node);
        void update_info();
        bool better_than(const Node *node); //TODO: this is a terrible name
        bool should_report(const AlnParams &params);

        size_t seed_len();
        size_t match_len();
        bool is_valid();

        
        void invalidate(std::vector<Node *> *old_nodes, bool delete_source);
        bool remove_child(Node *child);
        
        size_t add_child(Node *child);
        void print() const;
    };
    typedef std::pair<Node *, Node::Type> parent_ptr;


    public:

    const SdslFMI &fmi_;
    Range *kmer_ranges_;

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

    GraphAligner(const SdslFMI &fmi, 
              const AlnParams &aln_params,
              const std::string &label);

    ~GraphAligner();

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

