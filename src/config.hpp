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

#ifndef _INCL_CONF
#define _INCL_CONF

#include <iostream>
#include <vector>
#include <cfloat>
#include "dtw.hpp"
#include "mapper.hpp"
#include "fast5_reader.hpp"
#include "toplevel_prms.hpp"
#include "toml.hpp"

#ifdef PYBIND
#include "pybind11/pybind11.h"
#endif

#define INDEX_SUFF ".uncl"

const std::string ACTIVE_STRS[] = {"full", "even", "odd"};
const std::string MODE_STRS[] = {"deplete", "enrich"};

#define GET_SET(T, N) T get_##N() { return N; } \
                      void set_##N(const T &v) { N = v; }

#define GET_SET_EXTERN(T, E, N) T get_##N() { return E.N; } \
                                void set_##N(const T &v) { E.N = v; } 

#define GET_TOML_EXTERN(T, V, S) GET_NAMED_TOML(T, S.V, #V)
#define GET_TOML(T, V) GET_NAMED_TOML(T, V, #V)
#define GET_NAMED_TOML(T, V, N) { \
    if (subconf.contains(N)) \
        V = toml::find<T>(subconf,N); \
}

class Config {

    public:

    enum class Mode {MAP, REALTIME, SIM, MAP_ORD, UNDEF};

    Mode mode;
    u16 threads;
    bool cprof;

    Mapper::Params        mapper         = Mapper::PRMS;
    ReadBuffer::Params    read_buffer    = ReadBuffer::PRMS;
    Normalizer::Params    normalizer     = Normalizer::PRMS_DEF;
    EventDetector::Params event_detector = EventDetector::PRMS_DEF;
    EventProfiler::Params event_profiler = EventProfiler::PRMS_DEF;
    SeedTracker::Params   seed_tracker   = SeedTracker::PRMS_DEF;
    Fast5Reader::Params   fast5_reader   = Fast5Reader::PRMS_DEF;
    PoreModelParams       pore_model     = PORE_MODEL_PRMS_DEF;

    RealtimeParams        realtime       = REALTIME_PRMS_DEF;
    SimulatorParams       simulator     = SIM_PRMS_DEF;
    MapOrdParams          map_pool_ord   = MAP_ORD_PRMS_DEF;
    DtwParams             dtw            = DTW_PRMS_DEF;

    Config() : mode(Mode::UNDEF), threads(1), cprof(false) {}

    Config(Mode m) : Config() {
        mode = m;
    }

    void set_mode_map_ord() {
        mapper.chunk_timeout = FLT_MAX;
        mapper.evt_timeout = FLT_MAX;
    }

    void set_r94_rna() {
        read_buffer.sample_rate = 3012;
        read_buffer.bp_per_sec  = 70;
        read_buffer.chunk_time  = 1.0;
        //read_buffer.max_chunks  = 20;
        read_buffer.seq_fwd     = false;

        //TODO make mapper store pore_model group
        mapper.pore_model  = "r94_rna";
        pore_model.name    = "r94_rna";
        pore_model.reverse    = true;

        mapper.min_seed_prob = -3.0;
        mapper.max_paths     = 10000;

        event_detector.window_length1 = 7;
        event_detector.window_length2 = 12;
        event_detector.threshold1     = 2.8;
        event_detector.threshold2     = 18.0;
    }

    void export_static() {
        mapper.normalizer = normalizer;
        mapper.event_detector = event_detector;
        mapper.event_profiler = event_profiler;
        mapper.seed_tracker = seed_tracker;

        Mapper::PRMS = mapper;
        ReadBuffer::PRMS = read_buffer;
    }

    void load_toml(const std::string &fname) {
        const auto config = toml::parse(fname);

        if (config.contains("global")) {
            const auto subconf = toml::find(config, "global");
            GET_TOML(u16, threads);
        }

        if (config.contains("realtime")) {
            const auto subconf = toml::find(config, "realtime");
            GET_TOML_EXTERN(std::string, host, realtime);
            GET_TOML_EXTERN(u16, port, realtime);
            GET_TOML_EXTERN(float, duration, realtime);
            GET_TOML_EXTERN(u32, min_active_reads, realtime);
            GET_TOML_EXTERN(u32, max_active_reads, realtime);

            if (subconf.contains("realtime_mode")) {
                std::string mode_str = toml::find<std::string>(subconf, "realtime_mode");
                for (u8 i = 0; i != (u8) RealtimeParams::Mode::NUM; i++) {
                    if (mode_str == MODE_STRS[i]) {
                        realtime.realtime_mode = (RealtimeParams::Mode) i;
                        break;
                    }
                }
            }

            if (subconf.contains("active_chs")) {
                std::string mode_str = toml::find<std::string>(subconf, "active_chs");
                for (u8 i = 0; i != (u8) RealtimeParams::ActiveChs::NUM; i++) {
                    if (mode_str == ACTIVE_STRS[i]) {
                        realtime.active_chs = (RealtimeParams::ActiveChs) i;
                        break;
                    }
                }
            }
        }

        //TODO rename subsection
        if (config.contains("simulator")) {
            const auto subconf = toml::find(config, "simulator");
            GET_TOML_EXTERN(std::string, ctl_seqsum, simulator);
            GET_TOML_EXTERN(std::string, unc_seqsum, simulator);
            GET_TOML_EXTERN(std::string, unc_paf, simulator);
            GET_TOML_EXTERN(float, sim_speed, simulator);
            GET_TOML_EXTERN(float, scan_time, simulator);
            GET_TOML_EXTERN(float, scan_intv_time, simulator);
            GET_TOML_EXTERN(float, ej_time, simulator);
        }

        if (config.contains("map_pool_ord")) {
            const auto subconf = toml::find(config, "map_pool_ord");
            GET_TOML_EXTERN(u32, min_active_reads, map_pool_ord);
        }

        if (config.contains("fast5_reader")) {
            const auto subconf = toml::find(config, "fast5_reader");

            GET_TOML_EXTERN(u32, max_buffer, fast5_reader);
            GET_TOML_EXTERN(u32, max_reads, fast5_reader);
            GET_TOML_EXTERN(std::vector<std::string>, fast5_files, fast5_reader);
            GET_TOML_EXTERN(std::vector<std::string>, read_filter, fast5_reader);
            GET_TOML_EXTERN(bool, load_bc,  fast5_reader)
        }

        if (config.contains("read_buffer")) {
            const auto subconf = toml::find(config, "read_buffer");

            GET_TOML_EXTERN(u32, start_chunk,  read_buffer)
            GET_TOML_EXTERN(u32, max_chunks,   read_buffer);
            GET_TOML_EXTERN(u16, bp_per_sec,   read_buffer);
            GET_TOML_EXTERN(u16, sample_rate,  read_buffer);
            GET_TOML_EXTERN(float, chunk_time, read_buffer);
            GET_TOML_EXTERN(u16, num_channels, read_buffer);
            GET_TOML_EXTERN(bool, skip_notempl,  read_buffer)
        }

        if (config.contains("mapper")) {
            const auto subconf = toml::find(config, "mapper");

            GET_TOML_EXTERN(u32, seed_len, mapper);
            GET_TOML_EXTERN(u32, min_rep_len, mapper);
            GET_TOML_EXTERN(u32, max_rep_copy, mapper);
            GET_TOML_EXTERN(u32, max_paths, mapper);
            GET_TOML_EXTERN(u32, max_consec_stay, mapper);
            GET_TOML_EXTERN(u32, max_events, mapper);
            GET_TOML_EXTERN(float, max_stay_frac, mapper);
            GET_TOML_EXTERN(float, min_seed_prob, mapper);
            GET_TOML_EXTERN(std::string, bwa_prefix, mapper);
            GET_TOML_EXTERN(std::string, idx_preset, mapper);
            GET_TOML_EXTERN(std::string, pore_model, mapper);
            GET_TOML_EXTERN(u16, evt_batch_size, mapper);
            GET_TOML_EXTERN(float, evt_timeout, mapper);
            GET_TOML_EXTERN(float, chunk_timeout, mapper);

            #ifdef DEBUG_OUT
            GET_TOML_EXTERN(std::string, dbg_prefix, mapper);
            #endif
        }

        if (config.contains("seed_tracker")) {
            const auto subconf = toml::find(config, "seed_tracker");

            GET_TOML_EXTERN(float, min_mean_conf,seed_tracker);
            GET_TOML_EXTERN(float, min_top_conf, seed_tracker);
            GET_TOML_EXTERN(u32, min_map_len, seed_tracker);
        }

        if (config.contains("normalizer")) {
            const auto subconf = toml::find(config, "normalizer");
            GET_TOML_EXTERN(u32, len, normalizer);
            GET_TOML_EXTERN(float, tgt_mean, normalizer);
            GET_TOML_EXTERN(float, tgt_stdv, normalizer);
        }

        if (config.contains("event_detector")) {
            const auto subconf = toml::find(config, "event_detector");
            GET_TOML_EXTERN(float, min_mean, event_detector);
            GET_TOML_EXTERN(float, max_mean, event_detector);
            GET_TOML_EXTERN(float, threshold1, event_detector);
            GET_TOML_EXTERN(float, threshold2, event_detector);
            GET_TOML_EXTERN(float, peak_height, event_detector);
            GET_TOML_EXTERN(u32, window_length1, event_detector);
            GET_TOML_EXTERN(u32, window_length2, event_detector);
        }

        if (config.contains("event_profiler")) {
            const auto subconf = toml::find(config, "event_profiler");
            GET_TOML_EXTERN(u32, win_len, event_profiler);
            GET_TOML_EXTERN(float, win_stdv_min, event_profiler);
        }

    }

    Config(const std::string &toml_file) : Config() {
        load_toml(toml_file);
    }

    GET_SET(u16, threads)


    //TODO define get<type, param>, set<type, param>, doc<type, param>
    //use templates
    #define GET_SET_DOC(C, T, N) \
        T get_##N() { return C##_prms.N; } \
        void set_##N(const T &v) { C##_prms.N = v; } \
        static const char * doc_##N() { return C##_docs.N; }

    GET_SET_EXTERN(std::vector<std::string>, fast5_reader, fast5_files)
    GET_SET_EXTERN(std::vector<std::string>, fast5_reader, read_filter)
    GET_SET_EXTERN(u32, fast5_reader, max_reads)
    GET_SET_EXTERN(u32, fast5_reader, max_buffer)

    GET_SET_EXTERN(std::string, realtime, host)
    GET_SET_EXTERN(u16, realtime, port)
    GET_SET_EXTERN(float, realtime, duration)
    GET_SET_EXTERN(u32, realtime, max_active_reads)
    GET_SET_EXTERN(RealtimeParams::ActiveChs, realtime, active_chs)
    GET_SET_EXTERN(RealtimeParams::Mode, realtime, realtime_mode)

    GET_SET_EXTERN(u32, event_profiler, win_len)
    GET_SET_EXTERN(float, event_profiler, win_stdv_min)
    
    GET_SET_EXTERN(std::string, mapper, bwa_prefix)
    GET_SET_EXTERN(std::string, mapper, idx_preset)
    GET_SET_EXTERN(std::string, mapper, pore_model)
    GET_SET_EXTERN(u32, mapper, max_events)
    GET_SET_EXTERN(u32, mapper, seed_len);

    #ifdef DEBUG_OUT
    GET_SET_EXTERN(std::string, mapper, dbg_prefix)
    #endif

    GET_SET_EXTERN(u16,   read_buffer, num_channels);
    GET_SET_EXTERN(u32,   read_buffer, start_chunk)
    GET_SET_EXTERN(u32,   read_buffer, max_chunks)
    GET_SET_EXTERN(float, read_buffer, chunk_time);
    GET_SET_EXTERN(float, read_buffer, sample_rate);


    GET_SET_EXTERN(std::string, simulator, ctl_seqsum);
    GET_SET_EXTERN(std::string, simulator, unc_seqsum);
    GET_SET_EXTERN(std::string, simulator, unc_paf);
    GET_SET_EXTERN(float, simulator, sim_speed);
    GET_SET_EXTERN(float, simulator, scan_time);
    GET_SET_EXTERN(float, simulator, scan_intv_time);
    GET_SET_EXTERN(float, simulator, ej_time);

    GET_SET_EXTERN(u32, map_pool_ord, min_active_reads);

    EventDetector::Params &get_event_detector_prms() {
        return event_detector;
    }

    void set_event_detector(const EventDetector::Params &p) {
        event_detector = p;
    }

    #ifdef PYBIND
    #define DEF_METH(P, D) c.def(#P, &Config::P, D);
    #define DEFPRP(P) c.def_property(#P, &Config::get_##P, &Config::set_##P);
    #define DEFPRP_DOC(P) c.def_property( \
                #P, \
                &Config::get_##P, \
                &Config::set_##P, \
                Config::doc_##P());

        
    static std::vector<std::string> _PARAM_GROUPS, _GLOBAL_PARAMS;

    static void pybind_defs(pybind11::class_<Config> &c) {
        c.def(pybind11::init<const std::string &>())
         .def(pybind11::init())
         .def("load_toml", &Config::load_toml);

        DEF_METH(set_r94_rna, "Sets parameters for RNA sequencing")
        DEF_METH(set_mode_map_ord, "Sets parameters for map_ord")
        DEF_METH(export_static, "Sets static parameters for ReadBuffer and Mapper")

        #define CONF_ATTR(P, D, V) \
            c.def_readwrite(#P, &Config::P, D); \
            V.push_back(std::string(#P));
        #define CONF_GROUP(P, D) CONF_ATTR(P, D, _PARAM_GROUPS)
        #define CONF_GLOBAL(P, D) CONF_ATTR(P, D, _GLOBAL_PARAMS)

        CONF_GLOBAL(threads, "Number of threads to use")

        CONF_GROUP(mapper, "Mapper parameters")
        CONF_GROUP(read_buffer, "ReadBuffer parameters")

        CONF_GROUP(normalizer, "") 
        CONF_GROUP(event_detector, "")
        CONF_GROUP(event_profiler, "")
        CONF_GROUP(seed_tracker, "")
        CONF_GROUP(fast5_reader, "") 
        CONF_GROUP(pore_model, "") 
        CONF_GROUP(realtime, "") 
        CONF_GROUP(simulator, "")
        CONF_GROUP(map_pool_ord, "")
        CONF_GROUP(dtw, "")

        c.def_readonly_static("_PARAM_GROUPS", &Config::_PARAM_GROUPS);
        c.def_readonly_static("_GLOBAL_PARAMS", &Config::_GLOBAL_PARAMS);

        DEFPRP(bwa_prefix)
        DEFPRP(idx_preset)
        DEFPRP(max_events)
        DEFPRP(seed_len);
        DEFPRP(chunk_time)

        #ifdef DEBUG_OUT
        DEFPRP(dbg_prefix)
        #endif

        DEFPRP(fast5_files)//_DOC
        DEFPRP(read_filter) //_DOC
        DEFPRP(max_reads) //_DOC
        DEFPRP(max_buffer)//_DOC

        DEFPRP(host)
        DEFPRP(port)
        DEFPRP(duration)
        DEFPRP(max_active_reads)
        DEFPRP(active_chs)
        DEFPRP(realtime_mode)

        DEFPRP(num_channels)
        DEFPRP(max_chunks)
        DEFPRP(sample_rate)

        DEFPRP(min_active_reads)

        DEFPRP(ctl_seqsum)
        DEFPRP(unc_seqsum)
        DEFPRP(unc_paf)
        DEFPRP(sim_speed)
        DEFPRP(scan_time)
        DEFPRP(scan_intv_time)
        DEFPRP(ej_time)

        #define PY_CONF_MODE(P) m.value(#P, Config::Mode::P);

        
        pybind11::enum_<Config::Mode> m(c, "Mode");
        PY_CONF_MODE(MAP);
        PY_CONF_MODE(REALTIME);
        PY_CONF_MODE(SIM);
        PY_CONF_MODE(MAP_ORD);
        //m.export_values();
    }
    #endif
};

#endif
