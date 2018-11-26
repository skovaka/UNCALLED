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

#ifndef CHANNEL_POOL_HPP
#define CHANNEL_POOL_HPP

#include <thread>
#include <vector>
#include <deque>
#include "mapper.hpp"

class Fast5Read {
    public:
    Fast5Read();
    Fast5Read(const std::string &filename);

    void swap(Fast5Read &r);
    float next_sig();
    bool empty() const;
    
    std::vector<float> signal;
    std::string name;
    u16 channel;
    u32 i;
};

class ChannelPool {
    public:
    ChannelPool(MapperParams &params, u16 nchannels, u16 nthreads);
    
    void add_fast5s(const std::vector<std::string> &new_fast5s);
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

        std::vector<Mapper> &mappers_;

        bool running_;

        //Corrasponding inputs/output
        std::vector< Fast5Read > inputs_;
        std::deque< ReadLoc > outputs_;

        std::thread thread_;
        std::mutex in_mtx_, out_mtx_;
    };

    //List of mappers - one for each channel
    std::vector<Mapper> mappers_;

    //Which thread each mapper is currently assigned to
    //AKA the Mapper mappings
    std::deque<u16> open_channels_;
    std::vector<bool> channel_busy_;

    //Store threads in order of # active mappers
    std::vector<u16> thread_ids_;
    
    std::vector<MapperThread> threads_;

    std::deque<std::string> filenames_;
    std::vector<Fast5Read> read_buffer_;
    Fast5Read next_read_;
    

};

#endif
