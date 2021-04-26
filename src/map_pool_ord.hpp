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

#ifndef _INCL_MAP_POOL_ORD
#define _INCL_MAP_POOL_ORD

#include <thread>
#include <vector>
#include <deque>
#include <unordered_set>
#include "conf.hpp"
#include "realtime_pool.hpp"
#include "fast5_reader.hpp"

class MapPoolOrd : public RealtimePool {
    public:

    //MapOrdParams PRMS;

    MapPoolOrd(Conf &conf);

    void add_fast5(const std::string &fname);
    void add_read(const std::string &id);
    void load_fast5s();

    std::vector<MapResult> update();
    void stop();
    bool running();

    void end_read(u16 channel);
    Fast5Read &get_read(u16 channel);

    private:

    bool try_add_chunk(ReadBuffer &chunk);
    bool is_read_finished(const ReadBuffer &r);

    Fast5Iter fast5s_;
    //RealtimePool pool_;

    u32 active_tgt_;

    using ChQueue = std::deque<Fast5Read>;
    std::vector<ChQueue> channels_;
    std::vector<u32> chunk_idx_;

    #ifdef PYDEBUG
    Mapper::Meta get_meta(u16 channel) {
        return mappers_[channel-1].meta_;
    }
    std::vector<bool> meta_hold_;
    #endif


    bool channels_empty_;

    #ifdef PYBIND

    #define PY_MAP_ORD_METH(P) c.def(#P, &MapPoolOrd::P);
    #define PY_MAP_ORD_PRM(P) p.def_readwrite(#P, &MapOrdParams::P);

    public:

    static void pybind_defs(pybind11::class_<MapPoolOrd, RealtimePool> &c) {
        c.def(pybind11::init<Conf &>());

        PY_MAP_ORD_METH(add_fast5)
        PY_MAP_ORD_METH(add_read)
        PY_MAP_ORD_METH(load_fast5s)
        PY_MAP_ORD_METH(update)
        PY_MAP_ORD_METH(stop)
        PY_MAP_ORD_METH(running)
        PY_MAP_ORD_METH(get_read);
        PY_MAP_ORD_METH(end_read);

        pybind11::class_<MapOrdParams> p(c, "MapOrdParams");
        PY_MAP_ORD_PRM(min_active_reads);

        #ifdef PYDEBUG
        PY_MAP_ORD_METH(get_meta)
        #endif
    }

    #endif
};


#endif
