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

    class PathBuffer {
        public:

        static unsigned char PROB_WIN_LEN, TYPE_WIN_LEN;

        unsigned short 
            length_, 
            consec_stays_;

        //Prob window head/tail index, type window head index
        unsigned char prtl_, prhd_, prlen_, tyhd_, tytl_, tylen_;

        float *prob_sums_;
        EventType *event_types_;

        unsigned short all_type_counts_[EventType::NUM_TYPES];
        unsigned char win_type_counts_[EventType::NUM_TYPES];

        Kmer prev_kmer_;
        bool sa_checked_;

        //Source constructor
        PathBuffer(Kmer kmer, float prob);

        //Sibling constructor
        PathBuffer(PathBuffer *a, Kmer kmer, float prob, EventType type);


        //Copy constructor
        //PathBuffer(PathBuffer *a);

        //Creates invalid node
        //PathBuffer();

        ~PathBuffer();

        void init_source(Kmer kmer, float prob);
        void init_from_sibling(PathBuffer *a, Kmer kmer, float prob, EventType type);
        void init_from_parent(PathBuffer *a, Kmer kmer, float prob, EventType type);
        void make_child(Kmer kmer, float prob, EventType type);
        void update_consec_stays();
        bool better_than_parent(const PathBuffer *a, float prob);
        bool better_than_sibling(const PathBuffer *a, float prob);
        bool should_report(const Range &r, const AlnParams &params, bool has_children);

        size_t event_len();
        size_t match_len() const;
        float mean_prob() const;
        float next_mean_prob();
        float next_mean_prob(float next_prob) const;

        void print() const;
    };

    public:

    const FMI &fmi_;
    Range *kmer_ranges_;

    AlnParams params_;
    std::string label_;
    
    std::vector<PathBuffer *> inactive_paths_;
    std::map<Range, PathBuffer *> prev_paths_, next_paths_;

    #ifdef VERBOSE_TIME
    float child_map_time_, child_add_time_, child_rpl_time_;
    #endif
    
    unsigned int cur_event_;
    Event prev_event_;

    LeafAligner(const FMI &fmi, 
              const AlnParams &aln_params,
              const std::string &label);

    ~LeafAligner();

    void new_read(size_t read_len);
    void reset();
    std::vector<Result> add_event(float *kmer_probs, std::ostream &seeds_out, std::ostream &timing_out);

    void print_graph(bool verbose);

    bool add_child(Range &range, 
                   PathBuffer *prev_path,
                   Kmer kmer,
                   float prob,
                   EventType type,
                   bool prev_is_sibling);

    size_t add_sources(const Range &range, Kmer kmer, float prob);

    std::vector<Result> pop_seeds(std::ostream &out); //Big result gathering loop
                                     //Probably split into a few methods
};


#endif

