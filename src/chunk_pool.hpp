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

class ReadStream {
    public:
    ReadStream(u16 chunk_len, const std::string &name);

    void set_name(const std::string &name);
    bool empty() const;

    std::vector<float> chunk_;
    std::string name_;
};

class ChunkPool {
    public:
    ChunkPool(MapperParams &params, u16 nchannels, u16 nthreads);
    
    void new_read(u16 ch, const std::string &name);
    bool add_chunk(u16 ch, std::vector<float> &chunk);

    std::vector<std::string> update();
    bool all_finished();
    void stop_all();

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

    //List of mappers - one for each channel
    std::vector<Mapper> mappers_;
    std::vector<MapperThread> threads_;
    std::vector<ReadStream> read_buffers_;

    std::vector<u16> buffer_queue_, active_queue_, out_chs_;
    std::vector<bool> channel_active_;

    //Store threads in order of # active mappers
    
};

#endif
