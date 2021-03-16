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
#include <mutex>
#include "bwa_index.hpp"
#include "normalizer.hpp"
#include "event_detector.hpp"
#include "event_profiler.hpp"
#include "pore_model.hpp"
#include "models.inl"
#include "seed_tracker.hpp"
#include "read_buffer.hpp"
#include "paf.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#endif

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
        u16 evt_batch_size;
        float evt_timeout;
        float chunk_timeout;

        std::string bwa_prefix;
        std::string idx_preset;
        std::string pore_model;

        SeedTracker::Params seed_tracker;
        Normalizer::Params normalizer;
        EventDetector::Params event_detector;
        EventProfiler::Params event_profiler;

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
    void reset();
    void set_failed();

    Paf map_read();

    void skip_events(u32 n);
    bool add_chunk(ReadBuffer &chunk);

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
    Paf get_paf() const;
    void deactivate();

    static const u8 EVENT_MOVE = 1,
                    EVENT_STAY = 0;
    static const std::array<u8,2> EVENT_TYPES;
    static std::array<u32,EVENT_TYPES.size()> EVENT_ADDS;
    static u32 PATH_MASK, PATH_TAIL_MOVE;

    //std::vector<bool> Mapper::unpack_moves(u64 moves, u8 length);


    using moves_t = u32;
    using moves_u8_t = std::array<u8,4>;
    static constexpr u8 U8_BITS = 8;

	static std::vector<bool> unpack_moves(u64 moves, u8 length);

    //static moves_u8_t moves_to_u8(moves_t moves) {
    //    moves_u8_t ret{0};
    //    for (size_t i = 0; i < ret.size(); i++) {
    //        for (size_t j = 0; j < U8_BITS; j++) {
    //            size_t k = U8_BITS * i + j;
    //            ret[i] |= ((moves >> k) & 1);
    //        }
    //    }
    //    return ret;
    //}

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

        std::vector<bool> get_moves() const;

        float prob_head() const;

        void free_buffers();
        void print() const;

        Range fm_range_;
        u8 length_,
           consec_stays_;

        //u8 stay_count_, move_count_;
        //u8 path_type_counts_[EVENT_TYPES.size()];
        moves_t event_moves_;

        u16 total_moves_;

        u16 kmer_;

        float seed_prob_;
        float *prob_sums_;

        bool sa_checked_;


        #ifdef PYDEBUG
        static u32 count_;
        u32 evt_, id_, parent_;
        #endif

        static void reset_count() {
        #ifdef PYDEBUG
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
    EventProfiler evt_prof_;
    Normalizer norm_;
    SeedTracker seed_tracker_;
    ReadBuffer read_;
    Paf out_;

    //u16 channel_;
    //u32 read_num_;
    bool chunk_processed_, last_chunk_, reset_;//, processing_, adding_;
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

    void dbg_seeds_out(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end, 
        u64 sa_start, 
        u32 ref_len
    );

    void dbg_paths_out();
    void dbg_events_out();
    void dbg_conf_out();

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
    std::deque<AnnoEvent> dbg_events_;
    void dbg_events_open();
    #endif

    #ifdef DEBUG_CONFIDENCE
    std::ofstream conf_out_;
    bool confident_mapped_;
    #endif

    #ifdef PYDEBUG
    public:
    //stores name, start, end, strand, event, path_buf, cluster
    //using DbgSeed = std::tuple<std::string, u64, u64, bool, u32, u32, u32>;

    struct DbgSeed {
        i32 ref_id;
        i64 start, end;
        bool fwd;
        u32 event, path, cluster;
    };

    struct DbgPath {
        u32 event, id, parent_event, parent_id;
        u64 fm_start, fm_end;
        u16 kmer;
        u32 length, total_moves;
        float match_prob, seed_prob;
        moves_t moves_pac;
    };

    void dbg_add_seed(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end, 
        u64 sa_start, 
        u32 ref_len
    );

    struct Debug {
        std::vector<Event> events;
        std::vector<EvtProf> evt_profs;
        std::vector<float> proc_signal;
        std::vector<DbgSeed> seeds;
        std::vector<DbgPath> paths;
        u32 conf_evt, conf_clust;
    };
    Debug dbg_;
    #endif

    #ifdef PYBIND

    #define PY_MAPPER_METH(N,D) map.def(#N, &Mapper::N, D);
    #define PY_MAPPER_PRM(P) prm.def_readwrite(#P, &Mapper::Params::P);
    #define PY_PATHBUF_PRM(P) path.def_readonly(#P, &Mapper::PathBuffer::P);

    public:

    static void pybind_defs(pybind11::class_<Mapper> &map) {

        map.def(pybind11::init());

        #ifdef PYDEBUG
        #define PY_MAPPER_DBG(P) dbg.def_readonly(#P, &Mapper::Debug::P);

        #define PY_DBG_ARRAY(T, P) dbg.def_property_readonly(#P, [](Mapper::Debug &d) \
             -> py::array_t<T> {return py::array_t<T>(d.P.size(), d.P.data());});                                                                         

        pybind11::class_<Debug> dbg(map, "Debug");
        PY_MAPPER_DBG(events)
        PY_MAPPER_DBG(evt_profs)
        PY_MAPPER_DBG(proc_signal)
        //PY_MAPPER_DBG(seeds)
        //PY_DBG_ARRAY(DbgSeed, seeds);
        
        dbg.def_property_readonly("seeds", [](Mapper::Debug &d) -> pybind11::array_t<DbgSeed> {
             return pybind11::array_t<DbgSeed>(d.seeds.size(), d.seeds.data());
        });
        dbg.def_property_readonly("paths", [](Mapper::Debug &d) -> pybind11::array_t<DbgPath> {
             return pybind11::array_t<DbgPath>(d.paths.size(), d.paths.data());
        });

        PY_MAPPER_DBG(conf_evt)
        PY_MAPPER_DBG(conf_clust)

        //map.def_static("moves_to_u8", Mapper::moves_to_u8);
        map.def_static("unpack_moves", Mapper::unpack_moves);

        PYBIND11_NUMPY_DTYPE(DbgSeed, ref_id, start, end, fwd, event, path, cluster);
        PYBIND11_NUMPY_DTYPE(DbgPath, event, id, parent_event, parent_id, fm_start, fm_end, kmer, length, total_moves, match_prob, seed_prob, moves_pac);

        #endif

        //event, id, parent_evt, parent_id
        //fm_start, fm_len, kmer, full_len
        //match_prob, seed_prob
        //moves
        pybind11::class_<PathBuffer> path(map, "PathBuffer");
        path.def_readonly("event", &Mapper::PathBuffer::evt_);
        path.def_readonly("id", &Mapper::PathBuffer::id_);
        path.def_readonly("parent", &Mapper::PathBuffer::parent_);
        path.def_readonly("kmer", &Mapper::PathBuffer::kmer_);
        path.def_readonly("total_moves", &Mapper::PathBuffer::total_moves_);
        path.def_readonly("seed_prob", &Mapper::PathBuffer::seed_prob_);
        path.def_property_readonly("match_prob", &Mapper::PathBuffer::prob_head);
        path.def_property_readonly("moves", &Mapper::PathBuffer::get_moves);
        path.def_property_readonly("fm_start", 
                [](Mapper::PathBuffer &p) -> u64 {return p.fm_range_.start_;}
        );
        path.def_property_readonly("fm_end", 
                [](Mapper::PathBuffer &p) -> u64 {return p.fm_range_.end_;}
        );

        //path.def_readonly("fm_len", &Mapper::PathBuffer::parent_);

        //PY_MAP_METH(add_read, "Adds a read ID to the read filter");
        //PY_MAP_METH(load_read_list, "Loads a list of read IDs from a text file to add to the read filter");
        //PY_MAP_METH(pop_read, "");
        //PY_MAP_METH(buffer_size, "");
        //PY_MAP_METH(fill_buffer, "");
        //PY_MAP_METH(all_buffered, "");
        //PY_MAP_METH(empty, "");

        pybind11::class_<Params> prm(map, "Params");
        PY_MAPPER_PRM(seed_len)
        PY_MAPPER_PRM(min_rep_len)
        PY_MAPPER_PRM(max_rep_copy)
        PY_MAPPER_PRM(max_paths)
        PY_MAPPER_PRM(max_consec_stay)
        PY_MAPPER_PRM(max_events)
        PY_MAPPER_PRM(max_stay_frac)
        PY_MAPPER_PRM(min_seed_prob)
        PY_MAPPER_PRM(evt_batch_size)
        PY_MAPPER_PRM(evt_timeout)
        PY_MAPPER_PRM(chunk_timeout)
        PY_MAPPER_PRM(bwa_prefix)
        PY_MAPPER_PRM(idx_preset)
        PY_MAPPER_PRM(pore_model)
        PY_MAPPER_PRM(seed_tracker)
        PY_MAPPER_PRM(normalizer)
        PY_MAPPER_PRM(event_detector)
        PY_MAPPER_PRM(event_profiler)

        #ifdef DEBUG_OUT
        PY_MAPPER_PRM(dbg_prefix)
        #endif
    }

    #endif

};


#endif
