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
#include "read_buffer.hpp"
#include "fast5_pool.hpp"
#include "mapper.hpp"
#include "params.hpp"

Fast5Pool::Fast5Pool(const std::string &fast5_list_fname, const std::string &read_filter_fname) {

    if (!read_filter_fname.empty()) {
        std::ifstream filter_file(read_filter_fname);
        std::string read;
        while (getline(filter_file, read)) {
            filter_.insert(read);
        }
    }

    fast5_list_ = load_fast5s(fast5_list_fname, reads_, 4000, filter_);

    threads_ = std::vector<MapperThread>(PARAMS.threads);

    for (u32 i = 0; i < PARAMS.threads; i++) {
        threads_[i].next_read_.swap(reads_.front());
        reads_.pop_front();
        threads_[i].in_buffered_ = true;
        threads_[i].start();
    }

}

std::vector<Paf> Fast5Pool::update() {
    std::vector<Paf> ret;

    for (u32 i = 0; i < PARAMS.threads; i++) {
        if (threads_[i].out_buffered_) {
            ret.push_back(threads_[i].paf_out_);
            threads_[i].out_buffered_ = false;
        }

        if (!threads_[i].in_buffered_) {
            if (reads_.empty()) threads_[i].running_ = false;
            else {
                threads_[i].next_read_.swap(reads_.front());
                reads_.pop_front();
                threads_[i].in_buffered_ = true;
            }
        }
    }

    if (reads_.size() < PARAMS.threads) {
        load_fast5s(fast5_list_, reads_, 4000, filter_);
    }

    return ret;
}

bool Fast5Pool::all_finished() {
    if (!reads_.empty()) return false;
    for (u16 i = 0; i < PARAMS.threads; i++) {
        if (!threads_[i].finished_) return false;
    }
    return true;
}

void Fast5Pool::stop_all() {
    reads_.clear();
    for (auto &t : threads_) {
        t.running_ = false;
        t.mapper_.request_reset();
        t.thread_.join();
    }
}

u16 Fast5Pool::MapperThread::THREAD_COUNT = 0;

Fast5Pool::MapperThread::MapperThread()
    : tid_(THREAD_COUNT++),
      running_(true),
      in_buffered_(false),
      out_buffered_(false),
      finished_(false) {
    
}

Fast5Pool::MapperThread::MapperThread(MapperThread &&mt) 
    : tid_(mt.tid_),
      running_(mt.running_),                                             
      in_buffered_(mt.in_buffered_), 
      out_buffered_(mt.in_buffered_), 
      finished_(mt.finished_),
      mapper_(),
      thread_(std::move(mt.thread_)) {}

void Fast5Pool::MapperThread::start() {
    thread_ = std::thread(&Fast5Pool::MapperThread::run, this);
}

void Fast5Pool::MapperThread::run() {
    while (running_) {

        while (!in_buffered_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        mapper_.new_read(next_read_);
        in_buffered_ = false;

        Paf p = mapper_.map_read();

        while (out_buffered_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        paf_out_ = p;
        out_buffered_ = true;
    }

    finished_ = true;
}
