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


class MapperMeta : public Mapper {
    public:

    static const u8 EVENT_MOVE = 1,
                    EVENT_STAY = 0;
    static const std::array<u8,2> EVENT_TYPES;
    static std::array<u32,EVENT_TYPES.size()> EVENT_ADDS;
    static u32 PATH_MASK, PATH_TAIL_MOVE;

    using moves_t = u32;
    using moves_u8_t = std::array<u8,4>;
    static constexpr u8 U8_BITS = 8;

	static std::vector<bool> unpack_moves(u64 moves, u8 length);

    private:

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
    void meta_open_all();
    void meta_close_all();

    void meta_seeds_out(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end, 
        u64 sa_start, 
        u32 ref_len
    );

    void meta_paths_out();
    void meta_events_out();
    void meta_conf_out();

    //Debug helper functions and variables
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

    struct MetaSeed {
        i32 ref_id;
        i64 start, end;
        bool fwd;
        u32 event, path, cluster;
    };

    struct MetaPath {
        u32 event, id, parent;
        u64 fm_start, fm_length;
        u16 kmer;
        u32 length, total_moves;
        float norm_pdf, seed_prob;
        moves_t moves_pac;
    };

    void meta_add_seed(
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
        std::vector<MetaSeed> seeds;
        std::vector<MetaPath> paths;
        u32 conf_evt, conf_clust;
    };
    Debug meta_;
    #endif

    #ifdef PYBIND

    #define PY_MAPPER_METH(N,D) map.def(#N, &Mapper::N, D);
    #define PY_MAPPER_PRM(P) prm.def_readwrite(#P, &Mapper::Params::P);
    #define PY_PATHBUF_PRM(P) path.def_readonly(#P, &Mapper::PathBuffer::P);

    public:

    static void pybind_defs(pybind11::class_<Mapper> &map) {

        map.def(pybind11::init());

        #ifdef PYDEBUG
        #define PY_MAPPER_DBG(P) meta.def_readonly(#P, &Mapper::Debug::P);

        #define PY_DBG_ARRAY(T, P) meta.def_property_readonly(#P, [](Mapper::Debug &d) \
             -> py::array_t<T> {return py::array_t<T>(d.P.size(), d.P.data());});                                                                         

        pybind11::class_<Debug> meta(map, "Debug");
        PY_MAPPER_DBG(events)
        PY_MAPPER_DBG(evt_profs)
        PY_MAPPER_DBG(proc_signal)

        PY_MAPPER_DBG(conf_evt)
        PY_MAPPER_DBG(conf_clust)
        
        meta.def_property_readonly("seeds", [](Mapper::Debug &d) -> pybind11::array_t<MetaSeed> {
             return pybind11::array_t<MetaSeed>(d.seeds.size(), d.seeds.data());
        });

        meta.def_property_readonly("paths", [](Mapper::Debug &d) -> pybind11::array_t<MetaPath> {
             return pybind11::array_t<MetaPath>(d.paths.size(), d.paths.data());
        });

        PYBIND11_NUMPY_DTYPE(MetaSeed, ref_id, start, end, fwd, event, path, cluster);
        PYBIND11_NUMPY_DTYPE(MetaPath, event, id, parent, fm_start, fm_length, kmer, 
                             length, total_moves, norm_pdf, seed_prob, moves_pac);


        //map.def_static("moves_to_u8", Mapper::moves_to_u8);
        map.def_static("unpack_moves", Mapper::unpack_moves);


        #endif

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
    }

    #endif

};


#endif
