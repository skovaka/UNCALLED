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

#ifndef CONF_HPP
#define CONF_HPP

#include <iostream>
#include <vector>
#include "timer.hpp"
#include "mapper.hpp"
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "event_detector.hpp"
#include "mapper.hpp"
#include "seed_tracker.hpp"
#include "read_buffer.hpp"
#include "pybind11/pybind11.h"

#define INDEX_SUFF ".uncl"

typedef struct {
    std::string fast5_list;
    std::string fast5_filter;
    u32 max_reads, max_buffer;
} Fast5Params;

const std::string ACTIVE_STRS[] = {"full", "even", "odd"};
const std::string MODE_STRS[] = {"deplete", "enrich"};

typedef struct {
    enum class ActiveChs {FULL, EVEN, ODD, NUM};
    ActiveChs active_chs;

    enum class Mode {DEPLETE, ENRICH, NUM};
    Mode realtime_mode;

    std::string host;
    u16 port;
    float duration;
} RealtimeParams;

#define GET_SET(T, N) T get_##N() { return N; } \
                      void set_##N(const T &v) { N = v; }

#define GET_SET_EXTERN(T, E, N) T get_##N() { return E.N; } \
                                void set_##N(const T &v) { E.N = v; }

#define DEFPRP(P) c.def_property(#P, &Conf::get_##P, &Conf::set_##P);

class Conf {
    public:
    Conf();
    Conf(const std::string &conf_file);

    void load_conf(const std::string &fname);

    void load_index_params();

    GET_SET(u16, threads)
    GET_SET(u16, num_channels)
    GET_SET(std::string, bwa_prefix)
    GET_SET(std::string, kmer_model)
    GET_SET(std::string, index_preset)

    GET_SET_EXTERN(std::string, fast5_prms, fast5_list)
    GET_SET_EXTERN(std::string, fast5_prms, fast5_filter)
    GET_SET_EXTERN(u32, fast5_prms, max_reads)
    GET_SET_EXTERN(u32, fast5_prms, max_buffer)

    GET_SET_EXTERN(std::string, realtime_prms, host)
    GET_SET_EXTERN(u16, realtime_prms, port)
    GET_SET_EXTERN(float, realtime_prms, duration)
    GET_SET_EXTERN(RealtimeParams::ActiveChs, realtime_prms, active_chs)
    GET_SET_EXTERN(RealtimeParams::Mode, realtime_prms, realtime_mode)

    GET_SET_EXTERN(u32, Mapper::PRMS, max_events)
    GET_SET_EXTERN(u32, Mapper::PRMS, max_chunks)
    GET_SET_EXTERN(float, Mapper::PRMS, chunk_size)

    static void add_pybind_vars(pybind11::class_<Conf> &c) {
        DEFPRP(threads)
        DEFPRP(num_channels)
        DEFPRP(bwa_prefix)
        DEFPRP(kmer_model)
        DEFPRP(index_preset)

        DEFPRP(fast5_list)
        DEFPRP(fast5_filter)
        DEFPRP(max_reads)
        DEFPRP(max_buffer)

        DEFPRP(host)
        DEFPRP(port)
        DEFPRP(duration)
        DEFPRP(active_chs)
        DEFPRP(realtime_mode)

        DEFPRP(max_events)
        DEFPRP(max_chunks)
        DEFPRP(chunk_size)
    }

    u16 threads, num_channels;
    std::string bwa_prefix;
    std::string kmer_model;
    std::string index_preset;

    Fast5Params fast5_prms;
    RealtimeParams realtime_prms;
};

#endif
