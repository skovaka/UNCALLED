#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include <iostream>
#include <vector>
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "seed_tracker.hpp"
#include "timer.hpp"

//#define DEBUG_TIME
//#define DEBUG_SEEDS

class AlnParams {
    public:
    AlnParams(const KmerModel &model,
              const BwaFMI &fmi,
              const std::string &probfn_fname,
              u32 seed_len, 
              u32 min_rep_len, 
              u32 max_rep_copy, 
              u32 max_paths, 
              u32 max_consec_stay,
              u32 min_aln_len_,
              float max_stay_frac,
              float min_seed_prob, 
              float min_mean_conf,
              float min_top_conf);
    
    float get_prob_thresh(u64 fm_length);
    float get_source_prob();

    const KmerModel &model_;
    const BwaFMI &fmi_;

    u32 seed_len_, 
        min_rep_len_,
        max_rep_copy_,
        max_paths_,
        max_consec_stay_,
        min_aln_len_;

    float max_stay_frac_,
          min_seed_prob_,
          min_mean_conf_,
          min_top_conf_;

    std::vector<u64> evpr_lengths_;
    std::vector<float> evpr_threshes_;
    std::vector<Range> kmer_fmranges_;
};


class Aligner {
    public:

    Aligner(const AlnParams &aln_params);

    ~Aligner();

    ReadAln add_event(const Event &event
                      #ifdef DEBUG_TIME
                      ,std::ostream &time_out
                      #endif
                      #ifdef DEBUG_SEEDS
                      ,std::ostream &seeds_out
                      #endif
                      );

    void reset();

    private:

    enum EventType { MATCH, STAY, NUM_TYPES };
    static const u8 TYPE_BITS = 1;

    class PathBuffer {
        public:
        PathBuffer();
        PathBuffer(const PathBuffer &p);

        void make_source(Range &range, 
                         u16 kmer, 
                         float prob);

        void make_child(PathBuffer &p, 
                        Range &range, 
                        u16 kmer, 
                        float prob, 
                        EventType type);

        void invalidate();
        bool is_valid() const;
        bool is_seed_valid(const AlnParams &params, 
                           bool has_children) const;

        u8 type_head() const;
        u8 type_tail() const;
        u8 match_len() const;

        void free_buffers();
        void print() const;

        static u8 MAX_PATH_LEN, TYPE_MASK;
        static u64 TYPE_ADDS[EventType::NUM_TYPES];

        Range fm_range_;
        u16 length_,
            kmer_,
            consec_stays_;

        float seed_prob_;
        float *prob_sums_;

        u64 event_types_;
        u8 path_type_counts_[EventType::NUM_TYPES];

        bool sa_checked_;
    };

    friend bool operator< (const PathBuffer &p1, const PathBuffer &p2);

    private:

    void update_seeds(PathBuffer &p, 
                      std::vector<Seed> &seeds, 
                      bool has_children);

    AlnParams params_;
    SeedTracker seed_tracker_;
    Range *kmer_ranges_;
    
    std::vector<float> kmer_probs_;
    std::vector<PathBuffer> prev_paths_, next_paths_;
    std::vector<bool> sources_added_;
    u32 prev_size_,
        event_i_;

    #ifdef DEBUG_TIME
    std::ostream &time_out_;
    double loop1_time_, fmrs_time_, fmsa_time_, 
           sort_time_, loop2_time_, fullsource_time_;
    #endif
    #ifdef DEBUG_SEEDS
    std::ostream &seeds_out_;
    #endif
};


#endif
