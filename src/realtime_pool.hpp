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

#ifndef CHUNK_POOL_HPP
#define CHUNK_POOL_HPP

#include <thread>
#include <vector>
#include <deque>
#include "mapper.hpp"
#include "conf.hpp"
#include "sync_out.hpp"

using MapResult = std::tuple<u16, u32, Paf>;

class RealtimePool {
    public:
    RealtimePool(Conf &conf);
    
    void start_timer();
    bool add_chunk(Chunk &chunk);
    bool try_add_chunk(Chunk &chunk);
    void end_read(u16 ch, u32 number);
    bool is_read_finished(const ReadBuffer &r);

    std::vector<MapResult> update();
    bool all_finished();
    void stop_all(); //TODO: just name stop

    private:

    class MapperThread {
        public:
        MapperThread(std::vector<Mapper> &mappers);
        MapperThread(MapperThread &&mt);

        void start();
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
    };

    void buffer_chunk(Chunk &c);

    RealtimeParams PRMS;


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
