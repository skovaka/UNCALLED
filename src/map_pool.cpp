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

#include <thread>
#include <chrono>
#include "map_pool.hpp"

MapPool::MapPool(Conf &conf)
    : fast5s_(conf.fast5_prms) {

    Mapper::load_static(conf.bwa_prefix, conf.kmer_model, conf.index_preset);

    threads_ = std::vector<MapperThread>(conf.threads);

    //fast5s_.fill_buffer();

    for (u32 i = 0; i < threads_.size(); i++) {
        //if (!fast5s_.empty()) {
        //    ReadBuffer r = fast5s_.pop_read();
        //    threads_[i].next_read_.swap(r);
        //    threads_[i].in_buffered_ = true;
        //}
        threads_[i].start();
    }
}

std::vector<Paf> MapPool::update() {
    std::vector<Paf> ret;

    fast5s_.fill_buffer();

    for (u32 i = 0; i < threads_.size(); i++) {
        if (threads_[i].out_buffered_) {
            ret.push_back(threads_[i].paf_out_);
            threads_[i].out_buffered_ = false;
        }


        if (!threads_[i].in_buffered_) {
            if (fast5s_.empty()) { 
                threads_[i].finished_ = true;
            } else {
                ReadBuffer r = fast5s_.pop_read();
                threads_[i].next_read_.swap(r);
                threads_[i].in_buffered_ = true;
            }
        }
    }

    return ret;
}

void MapPool::add_fast5(const std::string &fast5_name) {
    fast5s_.add_fast5(fast5_name);
}


bool MapPool::running() {
    for (u16 i = 0; i < threads_.size(); i++) {
        if (threads_[i].running_) return true;
    }
    return false;
}

void MapPool::stop() {
    #ifdef FM_PROFILER
    FMProfiler prof_combined;
    #endif

    //reads_.clear();
    for (auto &t : threads_) {
        t.stopped_ = true;
        t.mapper_.request_reset();
        t.thread_.join();

        #ifdef FM_PROFILER
        prof_combined.combine(t.mapper_.fm_profiler_);
        #endif
    }

    #ifdef FM_PROFILER
    prof_combined.write("query_counts.bed");
    #endif
}

u16 MapPool::MapperThread::THREAD_COUNT = 0;

MapPool::MapperThread::MapperThread()
    : tid_(THREAD_COUNT++),
      running_(true),
      stopped_(false),
      finished_(false),
      in_buffered_(false),
      out_buffered_(false) {
    
}

MapPool::MapperThread::MapperThread(MapperThread &&mt) 
    : tid_(mt.tid_),
      running_(mt.running_),
      stopped_(mt.stopped_),                                             
      finished_(mt.finished_),                                             
      in_buffered_(mt.in_buffered_), 
      out_buffered_(mt.in_buffered_), 
      mapper_(),
      thread_(std::move(mt.thread_)) {}

void MapPool::MapperThread::start() {
    thread_ = std::thread(&MapPool::MapperThread::run, this);
}

void MapPool::MapperThread::run() {
    running_ = true;

    while (!(finished_ || stopped_)) {
        while (!in_buffered_ && !(stopped_ || finished_)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        if (finished_ || stopped_) break;

        mapper_.new_read(next_read_);
        in_buffered_ = false;

        Paf p = mapper_.map_read();

        while (out_buffered_ && !stopped_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        paf_out_ = p;
        out_buffered_ = !stopped_;
    }

    while (out_buffered_ && !stopped_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    running_ = false;
}
