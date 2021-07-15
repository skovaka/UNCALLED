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
#include "ref_index.hpp"
#include "normalizer.hpp"
#include "event_detector.hpp"
#include "event_profiler.hpp"
#include "pore_model.hpp"
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

constexpr KmerLen KLEN = KmerLen::k5;

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
        std::string pore_model; //TODO use PoreModel::params

        SeedTracker::Params seed_tracker;
        Normalizer::Params normalizer;
        EventDetector::Params event_detector;
        EventProfiler::Params event_profiler;

        #ifdef DEBUG_OUT
        std::string meta_prefix;
        #endif
    } Params;

    static Params PRMS;

    //TODO PRIVATIZE
    static RefIndex<KLEN> fmi;
    static PoreModel<KLEN> model;
    static std::vector<float> prob_threshes_;

    static void load_static();
    static inline i64 get_fm_bin(i64 fmlen);

    enum class State { INACTIVE, MAPPING, SUCCESS, FAILURE };

    Mapper();
    Mapper(const Mapper &m);

    ~Mapper();

    float get_prob_thresh(i64 fmlen) const;
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

    using moves_t = u32;
    using moves_u8_t = std::array<u8,4>;
    static constexpr u8 U8_BITS = 8;

	static std::vector<bool> unpack_moves(u64 moves, u8 length);

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


    //Meta output functions
    //All will be empty if no DEBUG_* macros are defined
    void meta_open_all();
    void meta_close_all();

    void meta_seeds_out(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end, 
        i64 sa_start, 
        u32 ref_len
    );

    void meta_paths_out();
    void meta_events_out();
    void meta_conf_out();

    //Meta helper functions and variables
    #ifdef DEBUG_OUT
    void meta_open(std::ofstream &out, const std::string &suffix);
    bool meta_opened_;
    #endif

    #ifdef DEBUG_SEEDS
    std::ofstream seeds_out_;
    void meta_seeds_open();
    #endif

    #ifdef DEBUG_PATHS
    std::ofstream paths_out_;
    void meta_paths_open();
    #endif

    #ifdef DEBUG_EVENTS
    std::ofstream events_out_;
    std::deque<AnnoEvent> meta_events_;
    void meta_events_open();
    #endif

    #ifdef DEBUG_CONFIDENCE
    std::ofstream conf_out_;
    bool confident_mapped_;
    #endif

    #ifdef PYDEBUG
    public:
    //stores name, start, end, strand, event, path_buf, cluster
    //using MetaSeed = std::tuple<std::string, u64, u64, bool, u32, u32, u32>;

    struct MetaEvent {
        u32 start, length;
        float mean, stdv;
        float prof_stdv, prof_mean;
        bool mask;
        float norm_sig;
    };

    struct MetaSeed {
        i32 ref_id;
        i64 start, end;
        bool fwd;
        u32 event, path, cluster;
    };

    struct MetaPath {
        u32 event, id, parent;
        i64 fm_start, fm_length;
        u16 kmer;
        u32 length, total_moves;
        float norm_pdf, seed_prob;
        moves_t moves_pac;
    };

    void meta_add_seed(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end, 
        i64 sa_start, 
        u32 ref_len
    );

    struct Meta {
        std::vector<MetaEvent> events;
        std::vector<u32> event_idxs;
        //std::vector<EvtProf> evt_profs;
        std::vector<float> proc_signal;
        std::vector<MetaSeed> seeds;
        std::vector<MetaPath> paths;
        u32 conf_evt, conf_clust;
    };
    Meta meta_;
    #endif

    #ifdef PYBIND

    #define PY_MAPPER_METH(N,D) map.def(#N, &Mapper::N, D);
    #define PY_MAPPER_PRM(P,D) prm.def_readwrite(#P, &Mapper::Params::P, D);
    #define PY_PATHBUF_PRM(P) path.def_readonly(#P, &Mapper::PathBuffer::P);

    public:

    static void pybind_defs(pybind11::class_<Mapper> &map) {

        map.def(pybind11::init());

        #ifdef PYDEBUG
        #define PY_MAPPER_DBG(P) meta.def_readonly(#P, &Mapper::Meta::P);

        #define PY_DBG_ARRAY(T, P) meta.def_property_readonly(#P, [](Mapper::Meta &d) \
             -> py::array_t<T> {return py::array_t<T>(d.P.size(), d.P.data());});                                                                         

        pybind11::class_<Meta> meta(map, "Meta");
        //PY_MAPPER_DBG(events)
        //PY_MAPPER_DBG(evt_profs)
        //PY_MAPPER_DBG(proc_signal)

        PY_MAPPER_DBG(conf_evt)
        PY_MAPPER_DBG(conf_clust)

        meta.def_property_readonly("masked_event_idxs", [](Mapper::Meta &d) -> pybind11::array_t<u32> {
             return pybind11::array_t<u32>(d.event_idxs.size(), d.event_idxs.data());
        });

        meta.def_property_readonly("events", [](Mapper::Meta &d) -> pybind11::array_t<MetaEvent> {
             return pybind11::array_t<MetaEvent>(d.events.size(), d.events.data());
        });
        
        meta.def_property_readonly("seeds", [](Mapper::Meta &d) -> pybind11::array_t<MetaSeed> {
             return pybind11::array_t<MetaSeed>(d.seeds.size(), d.seeds.data());
        });

        meta.def_property_readonly("paths", [](Mapper::Meta &d) -> pybind11::array_t<MetaPath> {
             return pybind11::array_t<MetaPath>(d.paths.size(), d.paths.data());
        });

        PYBIND11_NUMPY_DTYPE(MetaEvent, start, length, mean, stdv, prof_stdv, prof_mean, mask, norm_sig);

        PYBIND11_NUMPY_DTYPE(MetaSeed, ref_id, start, end, fwd, event, path, cluster);

        PYBIND11_NUMPY_DTYPE(MetaPath, event, id, parent, fm_start, fm_length, kmer, 
                             length, total_moves, norm_pdf, seed_prob, moves_pac);


        //map.def_static("moves_to_u8", Mapper::moves_to_u8);
        map.def_static("unpack_moves", Mapper::unpack_moves);


        #endif

        pybind11::class_<Params> prm(map, "Params");
        PY_MAPPER_PRM(seed_len, "Minimum length of seed to consider adding to tracker")
        PY_MAPPER_PRM(min_rep_len, "Minimum length of repetitve seeds")
        PY_MAPPER_PRM(max_rep_copy, "Maximum number of seeds to create from one path")
        PY_MAPPER_PRM(max_paths, "Maximum number of paths to store per event")
        PY_MAPPER_PRM(max_consec_stay, "Maximum number of consecutive stays")
        PY_MAPPER_PRM(max_events, "Will give up on a read after this many events have been processed")
        PY_MAPPER_PRM(max_stay_frac, "Maximum fraction of stay events per seed")
        PY_MAPPER_PRM(min_seed_prob, "Minimum average log event/k-mer match probability per seed")
        PY_MAPPER_PRM(evt_batch_size, "Number of events to process per thread iteration")
        PY_MAPPER_PRM(evt_timeout, "Maximum time (ms) to spend per event (averaged over event batch)")
        PY_MAPPER_PRM(chunk_timeout, "Maximum time (ms) to wait for a new chunk before read ends in realtime mode")
        PY_MAPPER_PRM(bwa_prefix, "BWA prefix to map to. Must be processed by \"uncalled index\".")
        PY_MAPPER_PRM(idx_preset, "Index parameters to use")
        PY_MAPPER_PRM(pore_model, "Filename or preset string for pore model")
        PY_MAPPER_PRM(seed_tracker, "Seed tracker parameters")
        PY_MAPPER_PRM(normalizer, "Normalizer parameters")
        PY_MAPPER_PRM(event_detector, "EventDetector parameters")
        PY_MAPPER_PRM(event_profiler, "EventProfiler parameters")
    }

    #endif

};


#endif
