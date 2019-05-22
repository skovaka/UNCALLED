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
    if (PARAMS.threads != 1) {
        std::cerr << "Multi-threaded mapping currently disabled. Using single thread\n";
    }

    if (!read_filter_fname.empty()) {
        std::ifstream filter_file(read_filter_fname);
        std::string read;
        while (getline(filter_file, read)) {
            filter_.insert(read);
        }
    }

    fast5_list_ = load_fast5s(fast5_list_fname, reads_, 4000, filter_);
    std::cerr << fast5_list_.size() << " " << reads_.size() << "\n";
    mappers_.emplace_back();
}

std::vector<Paf> Fast5Pool::update() {
    std::vector<Paf> ret;
    if (reads_.empty()) return ret;
    
    mappers_[0].new_read(reads_.front());
    ret.push_back(mappers_[0].map_read());
    reads_.pop_front();

    if (reads_.empty() && !fast5_list_.empty()) {
        load_fast5s(fast5_list_, reads_, 4000, filter_);
        std::cerr << fast5_list_.size() << " " << reads_.size() << "\n";
    }

    

    return ret;
}

bool Fast5Pool::all_finished() {
    return reads_.empty() && fast5_list_.empty();
}

void Fast5Pool::stop_all() {
}

/*
Fast5Pool::MapperThread::MapperThread(const UncalledOpts &opts)
    : running_(true),
      aligning_(false),
      mapper_(opts) {}

Fast5Pool::MapperThread::MapperThread(MapperThread &&mt) 
    : running_(mt.running_),                                             
      aligning_(mt.aligning_),                                           
      mapper_(mt.mapper_),
      thread_(std::move(mt.thread_)) {}

void Fast5Pool::MapperThread::start() {
    thread_ = std::thread(&Fast5Pool::MapperThread::run, this);
}

void Fast5Pool::MapperThread::run() {
    std::string fast5_id;
    std::vector<float> fast5_signal;
    ReadLoc loc;
    while (running_) {
        if (signals_in_.empty()) {
            aligning_ = false;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }

        aligning_ = true;
        
        in_mtx_.lock(); 
        signals_in_.front().swap(fast5_signal);
        ids_in_.front().swap(fast5_id);
        ids_in_.pop_front();
        signals_in_.pop_front();
        in_mtx_.unlock();

        mapper_.new_read(fast5_id);
        loc = mapper_.add_samples(fast5_signal);

        out_mtx_.lock(); 
        locs_out_.push_back(loc);
        out_mtx_.unlock();
    }
}
*/
