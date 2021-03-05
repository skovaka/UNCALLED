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

#ifndef _INCL_REALTIME_POOL
#define _INCL_REALTIME_POOL

#include <thread>
#include <vector>
#include <deque>
#include "mapper.hpp"
#include "conf.hpp"

//TODO create struct, auto cast to tuple?
using MapResult = std::tuple<u16, u32, Paf, bool>;

class RealtimePool {
    public:

    //static Params const PRMS_DEF;
    RealtimeParams PRMS;

    RealtimePool(Conf &conf);
    
    bool add_chunk(Chunk &chunk);
    bool try_add_chunk(Chunk &chunk);
    void end_read(u16 ch, u32 number);
    bool is_read_finished(const ReadBuffer &r);

    bool is_stopped() {return stopped_;}

    std::vector<MapResult> update();
    bool all_finished();
    void stop_all(); //TODO: just name stop

    u32 active_count() const; 

    #ifdef PYBIND

    #define PY_REALTIME_METH(P) c.def(#P, &RealtimePool::P);
    #define PY_REALTIME_PRM(P) p.def_readwrite(#P, &RealtimeParams::P);
    #define PY_REALTIME_MODE(P) m.value(#P, RealtimeParams::Mode::P);
    #define PY_REALTIME_ACTIVE(P) a.value(#P, RealtimeParams::ActiveChs::P);

    static void pybind_defs(pybind11::class_<RealtimePool> &c) {
        c.def(pybind11::init<Conf &>());
        PY_REALTIME_METH(add_chunk);
        PY_REALTIME_METH(try_add_chunk);
        PY_REALTIME_METH(update);
        PY_REALTIME_METH(all_finished);
        PY_REALTIME_METH(stop_all);

        pybind11::class_<RealtimeParams> p(c, "RealtimeParams");
        PY_REALTIME_PRM(host);
        PY_REALTIME_PRM(port);
        PY_REALTIME_PRM(duration);
        PY_REALTIME_PRM(max_active_reads);
        PY_REALTIME_PRM(active_chs);
        PY_REALTIME_PRM(realtime_mode);

        pybind11::enum_<RealtimeParams::Mode> m(c, "RealtimeMode");
        PY_REALTIME_MODE(DEPLETE);
        PY_REALTIME_MODE(ENRICH);
        m.export_values();

        pybind11::enum_<RealtimeParams::ActiveChs> a(c, "ActiveChs");
        PY_REALTIME_ACTIVE(FULL);
        PY_REALTIME_ACTIVE(EVEN);
        PY_REALTIME_ACTIVE(ODD);
        a.export_values();
    }

    #endif

    private:

    class MapperThread {
        public:
        MapperThread(std::vector<Mapper> &mappers);
        MapperThread(MapperThread &&mt);

        void start();
        void stop();

        void run();

        u16 read_count() const;

        static u16 num_threads;
        u16 tid_;

        std::vector<Mapper> &mappers_;

        bool running_;

        //Corrasponding inputs/output
        std::vector< u16 > in_chs_, in_tmp_, 
                           out_chs_, out_tmp_,
                           active_chs_;
        std::mutex in_mtx_, out_mtx_;

        std::thread thread_;

        float mtx_time_;
    };
    
    //TODO move to different module?
    enum class MapAction {
        ENRICH,  //Keep if mapped (=)
        DEPLETE, //Eject if mapped (*)
        FWD,     //Keep if mapped on + strand, eject otherwise (+)
        REV      //Keep if mapped on - strand, eject otherwise (-_
    };
    void set_targets(const RefTargets &tgts, MapAction action);

    void buffer_chunk(Chunk &c);
    bool should_keep(u32 rid, const Paf &paf);

    bool stopped_;

    u32 active_count_;

    std::vector<MapAction> tgt_actions_;

    //List of mappers - one for each channel
    std::vector<Mapper> mappers_;
    std::vector<MapperThread> threads_;
    std::vector<Chunk> chunk_buffer_;

    std::vector<u16> buffer_queue_, out_chs_, active_queue_;
    //std::deque<u16> ;
    //std::vector<u16> active_queue_;

    Timer time_;
    //Store threads in order of # active mappers
};

#endif
