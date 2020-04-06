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

#ifndef MAP_POOL_HPP
#define MAP_POOL_HPP

#include <thread>
#include <vector>
#include <deque>
#include <unordered_set>
#include "conf.hpp"

//typedef struct {
//    std::string fast5_list;
//    std::string filter;
//    u32 max_reads, max_buffer;
//} Fast5Params;


class Fast5Loader {
    public:
    //Fast5Loader(std::string fast5_list_fname, u32 max_buffer, std::string read_filter_fname="", u32 max_reads=UINT_MAX);
    Fast5Loader(const Fast5Params &p);

    void add_fast5(const std::string &fast5_name);

    ReadBuffer pop_read();
    u32 buffered_count();
    u32 buffer_reads();
    bool empty();

    Fast5Params PRMS;

    private:
    enum Format {MULTI, SINGLE, UNKNOWN};
    static const std::string FMT_RAW_PATHS[], FMT_CH_PATHS[];

    bool open_next();

    u32 max_buffer_, total_buffered_, max_reads_;

    std::deque<std::string> fast5_list_;
    std::unordered_set<std::string> filter_;

    hdf5_tools::File open_fast5_;
    Format open_fmt_;
    std::deque<std::string> read_paths_;

    std::deque<ReadBuffer> buffered_reads_;
};

class Fast5Pool {
    public:

    Fast5Pool(Conf &conf);

    std::vector<Paf> update();

    bool running();
    void add_fast5(const std::string &fname);
    void stop();

    private:
    Fast5Loader fast5s_;
    //std::deque<std::string> fast5_list_;
    //std::deque<ReadBuffer> reads_;
    //std::unordered_set<std::string> filter_;

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
