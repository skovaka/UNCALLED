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
#include "fast5_pool.hpp"


const std::string Fast5Loader::FMT_RAW_PATHS[] = {
    "/Raw",      //MULTI
    "/Raw/Reads" //SINGLE
};

const std::string Fast5Loader::FMT_CH_PATHS[] = {
    "/channel_id",                //MULTI
    "/UniqueGlobalKey/channel_id" //SINGLE
};

Fast5Loader::Fast5Loader(const Fast5Params &p) : PRMS(p) {
    total_buffered_ = 0;
    
    //std::ifstream list_file(PRMS.fast5_list);
    //std::string fast5_name;
    //while (getline(list_file, fast5_name)) {
    //    fast5_list_.push_back(fast5_name);
    //}
    
    if (!PRMS.fast5_filter.empty()) {
        std::ifstream filter_file(PRMS.fast5_filter);
        std::string read_name;

        while (getline(filter_file, read_name) && (PRMS.max_reads == 0 || filter_.size() < PRMS.max_reads)) {
            filter_.insert(read_name);
        }
    }
}

void Fast5Loader::add_fast5(const std::string &fast5_name) {
    fast5_list_.push_back(fast5_name);
}

bool Fast5Loader::empty() {
    return buffered_reads_.empty() && fast5_list_.empty();
}

bool Fast5Loader::open_next() {

    read_paths_.clear();
    if (open_fast5_.is_open()) open_fast5_.close();
    if (fast5_list_.empty()) return false;

    open_fast5_.open(fast5_list_.front());
    fast5_list_.pop_front();

    open_fmt_ = Format::UNKNOWN;
    for (const std::string &s : open_fast5_.list_group("/")) {
        if (s == "Raw") {
            open_fmt_ = Format::SINGLE;
            break;
        }
    }
    if (open_fmt_ == Format::UNKNOWN) open_fmt_ = Format::MULTI; //TODO: add support for old multi format


    switch (open_fmt_) {
    case Format::SINGLE:
        read_paths_.push_back("");
        return true;

    case Format::MULTI:
        for (const std::string &read : open_fast5_.list_group("/")) {
            std::string id = read.substr(read.find('_')+1);
            if (filter_.empty() || filter_.count(id) > 0) {
                read_paths_.push_back("/"+read);
            }
        }
        return true;
    default:
        return false;
    }

    return false; 
}

u32 Fast5Loader::buffer_reads() {
    u32 count = 0;

    //TODO: max total default to max int
    while (buffered_reads_.size() < PRMS.max_buffer && (PRMS.max_reads == 0 || total_buffered_ < PRMS.max_reads)) {

        while (read_paths_.empty()) {
            if(!open_next()) break;
        }
        if (read_paths_.empty()) break;

        std::string raw_path = read_paths_.front() + FMT_RAW_PATHS[open_fmt_],
                    ch_path =  read_paths_.front() + FMT_CH_PATHS[open_fmt_];
        read_paths_.pop_front();

        buffered_reads_.emplace_back(open_fast5_, raw_path, ch_path);


        count++;
        total_buffered_++;
    }


    if (PRMS.max_reads != 0 && total_buffered_ >= PRMS.max_reads) {
        read_paths_.clear();
        fast5_list_.clear();
    }

    return count;
}

ReadBuffer Fast5Loader::pop_read() {
    //TODO: swap to speed up?
    ReadBuffer r = buffered_reads_.front();
    buffered_reads_.pop_front();
    return r;
}

u32 Fast5Loader::buffered_count() {
    return buffered_reads_.size();
}

Fast5Pool::Fast5Pool(Conf &conf)
    : fast5s_(conf.fast5_prms) {

    conf.load_index_params();
    Mapper::model = KmerModel(conf.kmer_model, true);
    Mapper::fmi = BwaFMI(conf.bwa_prefix, Mapper::model);

    //fast5s_.buffer_reads();

    threads_ = std::vector<MapperThread>(conf.threads);

    for (u32 i = 0; i < threads_.size(); i++) {
        if (!fast5s_.empty()) {
            ReadBuffer r = fast5s_.pop_read();
            threads_[i].next_read_.swap(r);
            threads_[i].in_buffered_ = true;
        }
        threads_[i].start();
    }
}

std::vector<Paf> Fast5Pool::update() {
    std::vector<Paf> ret;

    fast5s_.buffer_reads();

    for (u32 i = 0; i < threads_.size(); i++) {
        if (threads_[i].out_buffered_) {
            ret.push_back(threads_[i].paf_out_);
            threads_[i].out_buffered_ = false;
        }


        if (!threads_[i].in_buffered_) {
            if (fast5s_.empty()) { 
                threads_[i].running_ = false;
            } else {
                ReadBuffer r = fast5s_.pop_read();
                threads_[i].next_read_.swap(r);
                threads_[i].in_buffered_ = true;
            }
        }
    }


    return ret;
}

void Fast5Pool::add_fast5(const std::string &fast5_name) {
    fast5s_.add_fast5(fast5_name);
}


bool Fast5Pool::all_finished() {
    if (!fast5s_.empty()) return false;
    for (u16 i = 0; i < threads_.size(); i++) {
        if (!threads_[i].finished_) return false;
    }
    return true;
}

void Fast5Pool::stop_all() {
    #ifdef FM_PROFILER
    FMProfiler prof_combined;
    #endif

    //reads_.clear();
    for (auto &t : threads_) {
	//t.running_ = false;
        t.stopped_ = true;
        //t.out_buffered_ = false;
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

u16 Fast5Pool::MapperThread::THREAD_COUNT = 0;

Fast5Pool::MapperThread::MapperThread()
    : tid_(THREAD_COUNT++),
      running_(true),
      stopped_(false),
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
    while (running_ && !stopped_) {
        while (!in_buffered_ && running_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        if (!running_) break;

        mapper_.new_read(next_read_);
        in_buffered_ = false;

        Paf p = mapper_.map_read();

        while (out_buffered_ && running_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        paf_out_ = p;
        out_buffered_ = !stopped_;
    }

    while (out_buffered_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    finished_ = true;
}
