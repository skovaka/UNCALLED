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

#define GET_TOML_EXTERN(C, T, V, S) GET_NAMED_TOML(C, T, S.V, #V)
#define GET_TOML(C, T, V) GET_NAMED_TOML(C, T, V, #V)
#define GET_NAMED_TOML(C, T, V, N) { \
    if (C.contains(N)) \
        V = toml::find<T>(C,N); \
}

class Conf {

    public:

    enum class Mode {MAP, REALTIME, SIM, MAP_ORD, UNDEF};

    Mode mode;
    u16 threads;

    Mapper::Params &mapper_prms = Mapper::PRMS;
    EventDetector::Params &event_prms = mapper_prms.event_prms;
    SeedTracker::Params &seed_prms = mapper_prms.seed_prms;

    ReadBuffer::Params &read_prms = ReadBuffer::PRMS;

    Fast5Reader::Params fast5_prms = Fast5Reader::PRMS_DEF;
    static constexpr Fast5Reader::Docstrs fast5_docs = Fast5Reader::DOCSTRS;

    RealtimeParams realtime_prms = REALTIME_PRMS_DEF;
    SimParams sim_prms = SIM_PRMS_DEF;
    MapOrdParams map_ord_prms = MAP_ORD_PRMS_DEF;

    Conf() : mode(Mode::UNDEF), threads(1) {}

    Conf(Mode m) : Conf() {
        mode = m;

        switch (mode) {

        case Mode::MAP_ORD:
            mapper_prms.chunk_timeout = FLT_MAX;
            mapper_prms.evt_timeout = FLT_MAX;
            break;

        default:
            break;
        }
    }

    void load_toml(const std::string &fname) {
        const auto conf = toml::parse(fname);

        if (conf.contains("global")) {
            const auto subconf = toml::find(conf, "global");
            GET_TOML(subconf, u16, threads);
        }

        if (conf.contains("realtime")) {
            const auto subconf = toml::find(conf, "realtime");
            GET_TOML_EXTERN(subconf, std::string, host, realtime_prms);
            GET_TOML_EXTERN(subconf, u16, port, realtime_prms);
            GET_TOML_EXTERN(subconf, float, duration, realtime_prms);
            GET_TOML_EXTERN(subconf, u32, max_active_reads, realtime_prms);

            if (subconf.contains("realtime_mode")) {
                std::string mode_str = toml::find<std::string>(subconf, "realtime_mode");
                for (u8 i = 0; i != (u8) RealtimeParams::Mode::NUM; i++) {
                    if (mode_str == MODE_STRS[i]) {
                        realtime_prms.realtime_mode = (RealtimeParams::Mode) i;
                        break;
                    }
                }
            }

            if (subconf.contains("active_chs")) {
                std::string mode_str = toml::find<std::string>(subconf, "active_chs");
                for (u8 i = 0; i != (u8) RealtimeParams::ActiveChs::NUM; i++) {
                    if (mode_str == ACTIVE_STRS[i]) {
                        realtime_prms.active_chs = (RealtimeParams::ActiveChs) i;
                        break;
                    }
                }
            }
        }

        //TODO rename subsection
        if (conf.contains("simulator")) {
            const auto subconf = toml::find(conf, "simulator");
            GET_TOML_EXTERN(subconf, std::string, ctl_seqsum, sim_prms);
            GET_TOML_EXTERN(subconf, std::string, unc_seqsum, sim_prms);
            GET_TOML_EXTERN(subconf, std::string, unc_paf, sim_prms);
            GET_TOML_EXTERN(subconf, float, sim_speed, sim_prms);
            GET_TOML_EXTERN(subconf, float, scan_time, sim_prms);
            GET_TOML_EXTERN(subconf, float, scan_intv_time, sim_prms);
            GET_TOML_EXTERN(subconf, float, ej_time, sim_prms);
            GET_TOML_EXTERN(subconf, u32, min_ch_reads, sim_prms);
        }

        if (conf.contains("map_ord")) {
            const auto subconf = toml::find(conf, "map_ord");
            GET_TOML_EXTERN(subconf, u32, min_active_reads, map_ord_prms);
        }

        if (conf.contains("fast5_reader")) {
            const auto subconf = toml::find(conf, "fast5_reader");

            GET_TOML_EXTERN(subconf, u32, max_buffer, fast5_prms);
            GET_TOML_EXTERN(subconf, u32, max_reads, fast5_prms);
            GET_TOML_EXTERN(subconf, std::string, fast5_list, fast5_prms);
            GET_TOML_EXTERN(subconf, std::string, read_list, fast5_prms);
        }

        if (conf.contains("reads")) {
            const auto subconf = toml::find(conf, "reads");

            GET_TOML_EXTERN(subconf, u32, max_chunks,   read_prms);
            GET_TOML_EXTERN(subconf, u16, bp_per_sec,   read_prms);
            GET_TOML_EXTERN(subconf, u16, sample_rate,  read_prms);
            GET_TOML_EXTERN(subconf, float, chunk_time, read_prms);
            GET_TOML_EXTERN(subconf, u16, num_channels, read_prms);
        }

        if (conf.contains("mapper")) {
            const auto subconf = toml::find(conf, "mapper");

            GET_TOML_EXTERN(subconf, u32, seed_len, mapper_prms);
            GET_TOML_EXTERN(subconf, u32, min_rep_len, mapper_prms);
            GET_TOML_EXTERN(subconf, u32, max_rep_copy, mapper_prms);
            GET_TOML_EXTERN(subconf, u32, max_paths, mapper_prms);
            GET_TOML_EXTERN(subconf, u32, max_consec_stay, mapper_prms);
            GET_TOML_EXTERN(subconf, u32, max_events, mapper_prms);
            GET_TOML_EXTERN(subconf, float, max_stay_frac, mapper_prms);
            GET_TOML_EXTERN(subconf, float, min_seed_prob, mapper_prms);
            GET_TOML_EXTERN(subconf, std::string, bwa_prefix, mapper_prms);
            GET_TOML_EXTERN(subconf, std::string, idx_preset, mapper_prms);
            GET_TOML_EXTERN(subconf, u32, evt_buffer_len, mapper_prms);
            GET_TOML_EXTERN(subconf, u16, evt_batch_size, mapper_prms);
            GET_TOML_EXTERN(subconf, float, evt_timeout, mapper_prms);
            GET_TOML_EXTERN(subconf, float, chunk_timeout, mapper_prms);

            #ifdef DEBUG_OUT
            GET_TOML_EXTERN(subconf, std::string, dbg_prefix, mapper_prms);
            #endif
        }

        if (conf.contains("seed_tracker")) {
            const auto subconf = toml::find(conf, "seed_tracker");

            GET_TOML_EXTERN(subconf, float, min_mean_conf,seed_prms);
            GET_TOML_EXTERN(subconf, float, min_top_conf, seed_prms);
            GET_TOML_EXTERN(subconf, u32, min_map_len, seed_prms);
        }

        if (conf.contains("event_detector")) {
            const auto subconf = toml::find(conf, "event_detector");
            GET_TOML_EXTERN(subconf, float, min_mean, event_prms);
            GET_TOML_EXTERN(subconf, float, max_mean, event_prms);
            GET_TOML_EXTERN(subconf, float, threshold1, event_prms);
            GET_TOML_EXTERN(subconf, float, threshold2, event_prms);
            GET_TOML_EXTERN(subconf, float, peak_height, event_prms);
            GET_TOML_EXTERN(subconf, u32, window_length1, event_prms);
            GET_TOML_EXTERN(subconf, u32, window_length2, event_prms);
        }

    }

    Conf(const std::string &toml_file) : Conf() {
        load_toml(toml_file);
    }

    GET_SET(u16, threads)


    //TODO define get<type, param>, set<type, param>, doc<type, param>
    //use templates
    #define GET_SET_DOC(C, T, N) \
        T get_##N() { return C##_prms.N; } \
        void set_##N(const T &v) { C##_prms.N = v; } \
        static const char * doc_##N() { return C##_docs.N; }

    GET_SET_DOC(fast5, std::string, fast5_list)
    GET_SET_DOC(fast5, std::string, read_list)
    GET_SET_DOC(fast5, u32, max_reads)
    GET_SET_DOC(fast5, u32, max_buffer)

    GET_SET_EXTERN(std::string, realtime_prms, host)
    GET_SET_EXTERN(u16, realtime_prms, port)
    GET_SET_EXTERN(float, realtime_prms, duration)
    GET_SET_EXTERN(u32, realtime_prms, max_active_reads)
    GET_SET_EXTERN(RealtimeParams::ActiveChs, realtime_prms, active_chs)
    GET_SET_EXTERN(RealtimeParams::Mode, realtime_prms, realtime_mode)

    GET_SET_EXTERN(std::string, mapper_prms, bwa_prefix)
    GET_SET_EXTERN(std::string, mapper_prms, idx_preset)
    GET_SET_EXTERN(u32, mapper_prms, max_events)

    #ifdef DEBUG_OUT
    GET_SET_EXTERN(std::string, mapper_prms, dbg_prefix)
    #endif

    GET_SET_EXTERN(u16,   read_prms, num_channels);
    GET_SET_EXTERN(u32,   read_prms, max_chunks)
    GET_SET_EXTERN(float, read_prms, chunk_time);
    GET_SET_EXTERN(float, read_prms, sample_rate);


    GET_SET_EXTERN(std::string, sim_prms, ctl_seqsum);
    GET_SET_EXTERN(std::string, sim_prms, unc_seqsum);
    GET_SET_EXTERN(std::string, sim_prms, unc_paf);
    GET_SET_EXTERN(float, sim_prms, sim_speed);
    GET_SET_EXTERN(float, sim_prms, scan_time);
    GET_SET_EXTERN(float, sim_prms, scan_intv_time);
    GET_SET_EXTERN(float, sim_prms, ej_time);
    GET_SET_EXTERN(u32, sim_prms, min_ch_reads);

    GET_SET_EXTERN(u32, map_ord_prms, min_active_reads);

    #ifdef PYBIND
    #define DEFPRP(P) c.def_property(#P, &Conf::get_##P, &Conf::set_##P);
    #define DEFPRP_DOC(P) c.def_property( \
                #P, \
                &Conf::get_##P, \
                &Conf::set_##P, \
                Conf::doc_##P());

    static void pybind_defs(pybind11::class_<Conf> &c) {
        c.def(pybind11::init<const std::string &>())
         .def(pybind11::init())
         .def("load_toml", &Conf::load_toml);

        DEFPRP(threads)

        DEFPRP(bwa_prefix)
        DEFPRP(idx_preset)
        DEFPRP(max_events)
        DEFPRP(chunk_time)

        #ifdef DEBUG_OUT
        DEFPRP(dbg_prefix)
        #endif

        DEFPRP_DOC(fast5_list)
        DEFPRP_DOC(read_list)
        DEFPRP_DOC(max_reads)
        DEFPRP_DOC(max_buffer)

        DEFPRP(host)
        DEFPRP(port)
        DEFPRP(duration)
        DEFPRP(max_active_reads)
        DEFPRP(active_chs)
        DEFPRP(realtime_mode)

        DEFPRP(num_channels)
        DEFPRP(max_chunks)

        DEFPRP(min_active_reads)

        DEFPRP(ctl_seqsum)
        DEFPRP(unc_seqsum)
        DEFPRP(unc_paf)
        DEFPRP(sim_speed)
        DEFPRP(scan_time)
        DEFPRP(scan_intv_time)
        DEFPRP(ej_time)
        DEFPRP(min_ch_reads)
    }
    #endif
};

#endif
