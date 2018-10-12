#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include <iostream>
#include <vector>
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"

//#define VERBOSE_TIME

//TODO: should be named seed
class Result {
    public:

    Result(unsigned int read_end, 
           unsigned int seed_len, 
           float prob, 
           unsigned int ref_start = 0, 
           unsigned int ref_end = 0);

    void set_ref_range(unsigned int end, unsigned int length);
    void print(std::ostream &out);

    Range read_range_, ref_range_;
    float seed_prob_;

    #ifdef DEBUG_PROB
    float min_evt_prob_;
    #endif
};

class AlnParams {
    public:
    const KmerModel &model_;

    unsigned int path_win_len_, 
                 min_rep_len_,
                 max_rep_copy_,
                 max_paths_;

    float max_stay_frac_;

    unsigned int max_consec_stay_, 
                 max_ignores_,
                 max_skips_;

    float window_prob_;

    std::vector<unsigned int> evpr_lengths_;
    std::vector<float> evpr_threshes_;

    AlnParams(const KmerModel &model,
              unsigned int path_win_len, 
              unsigned int min_rep_len, 
              unsigned int max_rep_copy, 
              unsigned int max_paths, 
              float max_stay_frac,
              unsigned int max_consec_stay,
              unsigned int max_ignores, 
              unsigned int max_skips,
              const std::string event_probs,
              float window_prob);           
    
    unsigned int nucl_to_events(unsigned int n);
    float get_prob_thresh(unsigned int fm_length);
    float get_source_prob();
};


class Aligner {
    enum EventType { MATCH, STAY, SKIP, IGNORE, NUM_TYPES };
    static const unsigned char TYPE_BITS = 2;

    class PathBuffer {
        public:

        static unsigned char MAX_WIN_LEN, TYPE_MASK;
        static unsigned long TYPE_ADDS[EventType::NUM_TYPES];

        unsigned short 
            length_, 
            consec_stays_;

        float win_prob_;

        float *prob_sums_;
        unsigned long event_types_;

        //unsigned short all_type_counts_[EventType::NUM_TYPES];
        unsigned char win_type_counts_[EventType::NUM_TYPES];

        Range fm_range_;
        u16 kmer_;
        bool sa_checked_;

        //Initial constructor
        PathBuffer();

        //Source constructor
        //PathBuffer(Range range, u16 kmer, float prob);

        //Child constructor
        //PathBuffer(PathBuffer &p, Range range, u16 kmer, float prob, EventType type);

        //Copy constructor
        PathBuffer(const PathBuffer &p);

        //~PathBuffer();

        //void init_from_sibling(PathBuffer *a, u16 kmer, float prob, EventType type);
        //void init_from_parent(PathBuffer *a, u16 kmer, float prob, EventType type);
        void make_source(Range &range, u16 kmer, float prob);
        void make_child(PathBuffer &p, Range &range, u16 kmer, float prob, EventType type);

        void invalidate();
        bool is_valid();

        unsigned char win_len() const;
        bool win_full() const;

        void update_consec_stays();
        //bool better_than_parent(const PathBuffer *a, float prob);
        bool better_than(const PathBuffer &a);
        unsigned char type_head() const;
        unsigned char type_tail() const;
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

    const BwaFMI &fmi_;
    Range *kmer_ranges_;

    AlnParams params_;
    
    //std::vector<PathBuffer *> inactive_paths_;
    std::vector<PathBuffer> prev_paths_, next_paths_;
    std::vector<bool> sources_added_;
    size_t prev_size_;

    #ifdef VERBOSE_TIME
    double loop1_time_, fmrs_time_, fmsa_time_, 
           sort_time_, loop2_time_, fullsource_time_;
    #endif
    
    unsigned int cur_event_;
    Event prev_event_;

    Aligner(const BwaFMI &fmi, 
                   const AlnParams &aln_params);

    ~Aligner();

    void check_alignments(PathBuffer &p, std::vector<Result> &results, bool has_children, std::ostream &seeds_out);

    void new_read(size_t read_len);
    void reset();

    //TODO: do probs inside
    std::vector<Result> add_event(std::vector<float> kmer_probs, std::ostream &seeds_out, std::ostream &timing_out);

    void print_graph(bool verbose);

    bool add_child(Range &range, 
                   PathBuffer &prev_path,
                   u16 kmer,
                   float prob,
                   EventType type,
                   bool prev_is_sibling);

    size_t add_sources(const Range &range, u16 kmer, float prob);
};


#endif
