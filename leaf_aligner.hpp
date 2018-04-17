#ifndef LEAF_ALIGNER_HPP
#define LEAF_ALIGNER_HPP

#include "fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"
#include "aligner.hpp"
#include <list>
#include <iostream>

class LeafAligner : public Aligner {

    enum EventType { MATCH, STAY, SKIP, IGNORE, NUM_TYPES };

    //TODO: think of better name
    //it's not an alignment (no FM info) - keeps track of events
    class Alignment {
        public:

        static unsigned char PROB_WIN_LEN, TYPE_WIN_LEN;

        unsigned short 
            length_, 
            consec_stays_;

        //Prob window head/tail index, type window head index
        unsigned char prtl_, prhd_, prlen_, tyhd_, tytl_, tylen_;

        double *prob_sums_;
        EventType *event_types_;

        unsigned short all_type_counts_[EventType::NUM_TYPES];
        unsigned char win_type_counts_[EventType::NUM_TYPES];

        Kmer prev_kmer_;

        //Source constructor
        Alignment(Kmer kmer, double prob);

        //Sibling constructor
        Alignment(Alignment *a, Kmer kmer, double prob, EventType type);


        //Copy constructor
        //Alignment(Alignment *a);

        //Creates invalid node
        //Alignment();

        ~Alignment();

        void init_source(Kmer kmer, double prob);
        void init_from_sibling(Alignment *a, Kmer kmer, double prob, EventType type);
        void init_from_parent(Alignment *a, Kmer kmer, double prob, EventType type);
        void make_child(Kmer kmer, double prob, EventType type);
        void update_consec_stays();
        bool better_than_parent(const Alignment *a, double prob);
        bool better_than_sibling(const Alignment *a, double prob);
        bool should_report(const AlnParams &params);

        size_t event_len();
        size_t match_len();
        double mean_prob() const;
        double next_mean_prob();
        double next_mean_prob(double next_prob) const;

        void print() const;
    };

    public:

    const FMI &fmi_;
    Range *kmer_ranges_;

    AlnParams params_;
    std::string label_;
    
    std::vector<Alignment *> inactive_alns_;
    std::map<Range, Alignment *> prev_alns_, next_alns_;

    Timer timer;
    
    unsigned int cur_event_;
    Event prev_event_;

    LeafAligner(const FMI &fmi, 
              const AlnParams &aln_params,
              const std::string &label);

    ~LeafAligner();

    void new_read(size_t read_len);
    void reset();
    std::vector<Result> add_event(double *kmer_probs, std::ostream &out);

    void print_graph(bool verbose);

    bool add_child(Range &range, 
                   Alignment *prev_aln,
                   Kmer kmer,
                   double prob,
                   EventType type,
                   bool prev_is_sibling);

    size_t add_sources(const Range &range, Kmer kmer, double prob);

    std::vector<Result> pop_seeds(std::ostream &out); //Big result gathering loop
                                     //Probably split into a few methods
};


#endif

