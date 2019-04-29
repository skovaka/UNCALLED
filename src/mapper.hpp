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

#ifndef MAPPER_HPP
#define MAPPER_HPP

#include <iostream>
#include <vector>
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "normalizer.hpp"
#include "seed_tracker.hpp"
#include "chunk.hpp"
#include "timer.hpp"
#include "read_buffer.hpp"

//#define DEBUG_TIME
//#define DEBUG_SEEDS

class Mapper {
    public:

    enum State { INACTIVE, MAPPING, SUCCESS, FAILURE };

    Mapper();
    Mapper(const Mapper &m);

    ~Mapper();


    void new_read(const std::string &name, u16 channel=0, u32 number=0);
    void new_read(Chunk &c);
    bool end_read(u32 number);

    Paf map_fast5(const std::string &fast5_name);
    //ReadLoc add_samples(const std::vector<float> &samples);
    //bool add_sample(float s);

    void skip_events(u32 n);
    bool swap_chunk(Chunk &chunk);

    u16 process_chunk();
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
        bool is_seed_valid(bool has_children) const;

        u8 type_head() const;
        u8 type_tail() const;
        u8 match_len() const;

        void free_buffers();
        void print() const;

        static u8 MAX_PATH_LEN, TYPE_MASK;
        static u32 TYPE_ADDS[EventType::NUM_TYPES];

        Range fm_range_;
        u8 length_,
            consec_stays_;
        u16 kmer_;

        float seed_prob_;
        float *prob_sums_;

        u32 event_types_;
        u8 path_type_counts_[EventType::NUM_TYPES];

        bool sa_checked_;
    };

    friend bool operator< (const PathBuffer &p1, const PathBuffer &p2);

    private:

    bool add_event(float event);

    void update_seeds(PathBuffer &p, bool has_children);

    void set_ref_loc(const SeedGroup &seeds);

    EventDetector event_detector_;
    Normalizer norm_;
    SeedTracker seed_tracker_;
    ReadBuffer read_;

    //u16 channel_;
    //u32 read_num_;
    bool last_chunk_, reset_;
    State state_;
    std::vector<float> kmer_probs_;
    std::vector<PathBuffer> prev_paths_, next_paths_;
    std::vector<bool> sources_added_;
    u32 prev_size_,
        event_i_,
        chunk_i_;
    Timer timer_;

    #ifdef DEBUG_SEEDS
    std::ostream &seeds_out_;
    #endif
};


#endif
