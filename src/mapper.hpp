#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include <iostream>
#include <vector>
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "event_detector.hpp"
#include "seed_tracker.hpp"
#include "timer.hpp"

//#define DEBUG_TIME
//#define DEBUG_SEEDS

class MapperParams {
    public:
    MapperParams(const std::string &bwa_prefix,
                 const std::string &model_fname,
                 const std::string &probfn_fname,
                 u32 seed_len, 
                 u32 min_aln_len,
                 u32 min_rep_len, 
                 u32 max_rep_copy, 
                 u32 max_consec_stay,
                 u32 max_paths, 
                 u32 max_events_proc,
                 u32 evt_winlen1,
                 u32 evt_winlen2,
                 float evt_thresh1,
                 float evt_thresh2,
                 float evt_peak_height,
                 float evt_min_mean,
                 float evt_max_mean,
                 float max_stay_frac,
                 float min_seed_prob, 
                 float min_mean_conf,
                 float min_top_conf);
    
    float get_prob_thresh(u64 fm_length) const;
    float get_source_prob() const;

    BwaFMI fmi_;
    KmerModel model_;
    EventParams event_params_;

    u32 seed_len_, 
        min_rep_len_,
        max_rep_copy_,
        max_paths_,
        max_consec_stay_,
        min_aln_len_,
        max_events_proc_;

    float max_stay_frac_,
          min_seed_prob_,
          min_mean_conf_,
          min_top_conf_;

    std::vector<u64> evpr_lengths_;
    std::vector<float> evpr_threshes_;
    std::vector<Range> kmer_fmranges_;
};

class ReadLoc {
    public:
    ReadLoc();
    ReadLoc(const std::string &rd_name);

    bool set_ref_loc(const MapperParams &params, const SeedGroup &seeds);
    void set_read_len(u32 len);
    std::string str() const;
    bool is_valid() const; 
    friend std::ostream &operator<< (std::ostream &out, const ReadLoc &l);

    private:
    std::string rd_name_, rf_name_;
    u64 rd_st_, rd_en_, rd_len_,
        rf_st_, rf_en_, rf_len_;
    u16 match_count_;
    bool fwd_;
};

std::ostream &operator<< (std::ostream &out, const ReadLoc &l);

class Mapper {
    public:

    Mapper(MapperParams &map_params);

    ~Mapper();

    void new_read(const std::string &name);

    std::string map_fast5(const std::string &fast5_name);
    ReadLoc add_samples(const std::vector<float> &samples);
    //SeedGroup add_sample(float s);

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
        bool is_seed_valid(const MapperParams &params, 
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

    bool add_event(const Event &event
                      #ifdef DEBUG_TIME
                      ,std::ostream &time_out
                      #endif
                      #ifdef DEBUG_SEEDS
                      ,std::ostream &seeds_out
                      #endif
                      );

    void update_seeds(PathBuffer &p, 
                      std::vector<SeedGroup> &seeds, 
                      bool has_children);

    const MapperParams &params_;
    const KmerModel &model_;
    const BwaFMI &fmi_;
    EventDetector event_detector_;
    SeedTracker seed_tracker_;

    ReadLoc read_loc_;
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
