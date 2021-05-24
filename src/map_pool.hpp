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

#ifndef _INCL_MAP_POOL
#define _INCL_MAP_POOL

#include <thread>
#include <vector>
#include <deque>
#include <unordered_set>
#include "config.hpp"

class MapPool {
    public:

    MapPool(Config &config);

    std::vector<Paf> update();

    bool running();
    void add_fast5(const std::string &fname);
    void stop();

    #ifdef PYBIND
    #define PY_MAP_POOL_METH(P) c.def(#P, &MapPool::P);

    static void pybind_defs(pybind11::class_<MapPool> &c) {
        c.def(pybind11::init<Config &>());
        PY_MAP_POOL_METH(update);
        PY_MAP_POOL_METH(running);
        PY_MAP_POOL_METH(add_fast5);
        PY_MAP_POOL_METH(stop);
    }

    #endif

    private:
    Fast5Iter fast5s_;

    class MapperThread {
        public:
        MapperThread();
        MapperThread(MapperThread &&mt);

        void start();
        void run();

        static u16 THREAD_COUNT;

        u16 tid_;

        //running: run method has not ended
        //stopped: force stopped, like by keyboard interrupt
        //finished: no more reads left to process
        bool running_, stopped_, finished_, 
             in_buffered_, out_buffered_;

        Mapper mapper_;
        std::thread thread_;

        ReadBuffer next_read_;
        Paf paf_out_;

        std::mutex in_mtx_, out_mtx_;
    };

    std::vector<MapperThread> threads_;

};


#endif
