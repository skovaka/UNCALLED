#ifndef LEAF_ARR_ALIGNER_HPP
#define LEAF_ARR_ALIGNER_HPP

#include "fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"
#include "aligner.hpp"
#include <list>
#include <iostream>

class LeafArrAligner : public Aligner {

    enum EventType { MATCH, STAY, SKIP, IGNORE, NUM_TYPES };

    class PathBuffer {
        public:

        static unsigned char PROB_WIN_LEN, TYPE_WIN_LEN;

        unsigned short 
            length_, 
            consec_stays_;

        //Prob window head/tail index, type window head index
        unsigned char prhd_, prtl_, prlen_, tyhd_, tytl_, tylen_;
        float win_prob_;

        float *prob_sums_;
        EventType *event_types_;

        char space_filler[54];
        unsigned short all_type_counts_[EventType::NUM_TYPES];
        unsigned char win_type_counts_[EventType::NUM_TYPES];

        Range fm_range_;
        Kmer kmer_;
        bool sa_checked_;

        //Initial constructor
        PathBuffer();

        //Source constructor
        //PathBuffer(Range range, Kmer kmer, float prob);

        //Child constructor
        //PathBuffer(PathBuffer &p, Range range, Kmer kmer, float prob, EventType type);

        //Copy constructor
        PathBuffer(const PathBuffer &p);

        //~PathBuffer();

        //void init_from_sibling(PathBuffer *a, Kmer kmer, float prob, EventType type);
        //void init_from_parent(PathBuffer *a, Kmer kmer, float prob, EventType type);
        void make_source(Range &range, Kmer kmer, float prob);
        void make_child(PathBuffer &p, Range &range, Kmer kmer, float prob, EventType type);

        void invalidate();
        bool is_valid();

        void update_consec_stays();
        //bool better_than_parent(const PathBuffer *a, float prob);
        bool better_than(const PathBuffer &a);
        bool should_report(const AlnParams &params, bool has_children);

        size_t event_len();
        size_t match_len() const;
        float mean_prob() const;
        float next_mean_prob();
        float next_mean_prob(float next_prob) const;

        void replace(const PathBuffer &r);
        //PathBuffer& operator=(const PathBuffer &r) = delete;

        void free_buffers();

        void print() const;
    };

    friend bool operator< (const PathBuffer &p1, const PathBuffer &p2);

    public:

    const FMI &fmi_;
    Range *kmer_ranges_;

    AlnParams params_;
    std::string label_;
    
    //std::vector<PathBuffer *> inactive_paths_;
    std::vector<PathBuffer> prev_paths_, next_paths_;
    std::vector<bool> sources_added_;
    size_t prev_size_;

    #ifdef VERBOSE_TIME
    float child_map_time_, child_add_time_, child_rpl_time_;
    #endif
    
    unsigned int cur_event_;
    Event prev_event_;

    LeafArrAligner(const FMI &fmi, 
              const AlnParams &aln_params,
              const std::string &label);

    ~LeafArrAligner();

    void check_alignments(PathBuffer &p, std::vector<Result> &results, bool has_children, std::ostream &seeds_out);

    void new_read(size_t read_len);
    void reset();
    std::vector<Result> add_event(float *kmer_probs, std::ostream &seeds_out, std::ostream &timing_out);

    void print_graph(bool verbose);

    bool add_child(Range &range, 
                   PathBuffer &prev_path,
                   Kmer kmer,
                   float prob,
                   EventType type,
                   bool prev_is_sibling);

    size_t add_sources(const Range &range, Kmer kmer, float prob);

    std::vector<Result> pop_seeds(std::ostream &out); //Big result gathering loop
                                     //Probably split into a few methods
};


#endif

