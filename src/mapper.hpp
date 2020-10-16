/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _INCL_MAPPER
#define _INCL_MAPPER

#include <iostream>
#include <vector>
#include "bwa_index.hpp"
#include "normalizer.hpp"
#include "event_detector.hpp"
#include "pore_model.hpp"
#include "seed_tracker.hpp"
#include "read_buffer.hpp"

const KmerLen KLEN = KmerLen::k5;

//#define DEBUG_TIME
//#define DEBUG_SEEDS

//TODO define as constant somewhere
//rematch "params" python module
#define INDEX_SUFF ".uncl"

class Mapper {
    public:

    typedef struct {
        //standard mapping
        u32 seed_len;
        u32 min_rep_len;
        u32 max_rep_copy;
        u32 max_paths;
        u32 max_consec_stay;
        u32 max_events;
        float max_stay_frac;
        float min_seed_prob;

        //realtime only
        u32 evt_buffer_len;
        u16 evt_batch_size;
        float evt_timeout;
        float chunk_timeout;

        std::string bwa_prefix, idx_preset;

        SeedTracker::Params seed_prms;
        EventDetector::Params event_prms;

        #ifdef DEBUG_OUT
        std::string dbg_prefix;
        #endif

    } Params;

    static Params PRMS;

    //TODO PRIVATIZE
    static BwaIndex<KLEN> fmi;
    static PoreModel<KLEN> model;
    static std::vector<float> prob_threshes_;

    static void load_static();
    static inline u64 get_fm_bin(u64 fmlen);

    enum class State { INACTIVE, MAPPING, SUCCESS, FAILURE };

    Mapper();
    Mapper(const Mapper &m);

    ~Mapper();

    float get_prob_thresh(u64 fmlen) const;
    float get_source_prob() const;
    u16 get_max_events() const;

    void new_read(ReadBuffer &r);
    void new_read(Chunk &c);
    void reset();
    void set_failed();

    Paf map_read();

    void skip_events(u32 n);
    bool add_chunk(Chunk &chunk);

    u32 event_to_bp(u32 evt_i, bool last=false) const;

    u32 events_mapped() const {return event_i_;}

    u16 process_chunk();
    bool chunk_mapped();
    bool map_chunk();
    bool is_chunk_processed() const;
    void request_reset();
    void end_reset();
    bool is_resetting();
    State get_state() const;

    u32 prev_unfinished(u32 next_number) const;

    bool finished() const;
    ReadBuffer &get_read();
    void deactivate();

    private:

    static const u8 EVENT_MOVE = 1,
                    EVENT_STAY = 0;
    static const std::array<u8,2> EVENT_TYPES;
    static std::array<u32,EVENT_TYPES.size()> EVENT_ADDS;
    static u32 PATH_MASK, PATH_TAIL_MOVE;
    //static u32 PATH_MASK;TODO popcount instead of store?

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
                        u8 event_type);

        void invalidate();
        bool is_valid() const;
        bool is_seed_valid(bool has_children) const;

        u8 type_head() const;
        u8 type_tail() const;
        u8 move_count() const;
        u8 stay_count() const;

        void free_buffers();
        void print() const;

        Range fm_range_;
        u8 length_,
           consec_stays_;

        //u8 stay_count_, move_count_;
        //u8 path_type_counts_[EVENT_TYPES.size()];
        u32 event_moves_;

        u16 total_move_len_;

        u16 kmer_;

        float seed_prob_;
        float *prob_sums_;

        bool sa_checked_;

        #ifdef DEBUG_OUT
        static u32 count_;
        u32 id_, parent_;
        #endif

        static void reset_count() {
        #ifdef DEBUG_OUT
            count_ = 0;
        #endif
        }
    };

    friend bool operator< (const PathBuffer &p1, const PathBuffer &p2);

    private:

    bool map_next();

    void update_seeds(PathBuffer &p, bool has_children);

    void set_ref_loc(const SeedCluster &seeds);


    EventDetector evdt_;
    Normalizer norm_;
    SeedTracker seed_tracker_;
    ReadBuffer read_;

    //u16 channel_;
    //u32 read_num_;
    bool last_chunk_, reset_;//, processing_, adding_;
    State state_;
    std::vector<float> kmer_probs_;
    std::vector<PathBuffer> prev_paths_, next_paths_;
    std::vector<bool> sources_added_;
    u32 prev_size_,
        event_i_,
        chunk_i_;
    Timer chunk_timer_, map_timer_;
    float map_time_, wait_time_;

    std::mutex chunk_mtx_;


    //Debug output functions
    //All will be empty if no DEBUG_* macros are defined
    void dbg_open_all();
    void dbg_close_all();
    void dbg_seeds_out(const PathBuffer &path, u32 clust, u64 ref_end, u32 evt_end);
    void dbg_paths_out();
    void dbg_events_out();

    //Debug helper functions and variables
    #ifdef DEBUG_OUT
    void dbg_open(std::ofstream &out, const std::string &suffix);
    bool dbg_opened_;
    #endif

    #ifdef DEBUG_SEEDS
    std::ofstream seeds_out_;
    void dbg_seeds_open();
    #endif

    #ifdef DEBUG_PATHS
    std::ofstream paths_out_;
    void dbg_paths_open();
    #endif

    #ifdef DEBUG_EVENTS
    std::ofstream events_out_;
    std::deque<Event> events_;
    void dbg_events_open();
    #endif

    #ifdef DEBUG_CONFIDENCE
    bool confident_mapped_;
    #endif
};


#endif
