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
#include <stdlib.h>
#include <time.h>
#include "realtime_pool.hpp"
#include "mapper.hpp"
#include "sync_out.hpp"

//SyncOut debug_out_(std::cout);

RealtimePool::RealtimePool(Conf &conf) :
    PRMS(conf.realtime_prms) {

    conf.load_index_params();
    Mapper::model = PoreModel<KLEN>(conf.kmer_model, true);
    Mapper::fmi.load_index(conf.bwa_prefix);

    for (u16 t = 0; t < conf.threads; t++) {
        threads_.emplace_back(mappers_);
    }

    mappers_.resize(conf.num_channels);
    chunk_buffer_.resize(conf.num_channels);
    buffer_queue_.reserve(conf.num_channels);
    active_queue_.reserve(conf.num_channels);


    for (u16 t = 0; t < conf.threads; t++) {
        threads_[t].start();
    }

    srand(time(NULL));
}

void RealtimePool::buffer_chunk(Chunk &c) {
    u16 ch = c.get_channel_idx();
    if (chunk_buffer_[ch].empty()) {
        buffer_queue_.push_back(ch);
    } else {
        //TODO: handle backlog (probably reset paths?)
        chunk_buffer_[ch].clear();
    }
    chunk_buffer_[ch].swap(c);
}


//Add chunk to master buffer
bool RealtimePool::add_chunk(Chunk &c) {
    u16 ch = c.get_channel_idx();

    //Check if previous read is still aligning
    //If so, tell thread to reset, store chunk in pool buffer
    if (mappers_[ch].prev_unfinished(c.get_number())) {
        mappers_[ch].request_reset();
        buffer_chunk(c);
        return true;

    //Previous alignment finished but mapper hasn't reset
    //Happens if update hasn't been called yet
    } else if (mappers_[ch].finished()) {
        if (mappers_[ch].get_read().number_ != c.get_number()){ 
            buffer_chunk(c);
        }
        return true;

    //Mapper inactive - need to reset graph and assign to thread
    } else if (mappers_[ch].get_state() == Mapper::State::INACTIVE) {
        mappers_[ch].new_read(c);
        active_queue_.push_back(ch);
        return true;

    }
    
    return mappers_[ch].add_chunk(c);
    

    return false;
}

bool RealtimePool::is_read_finished(const ReadBuffer &r) {
    u16 ch = r.get_channel_idx();
    return (mappers_[ch].finished() && 
            mappers_[ch].get_read().get_number() == r.get_number());
}

bool RealtimePool::try_add_chunk(Chunk &c) {
    u16 ch = c.get_channel_idx();

    //Chunk is empty if all read chunks were output
    if (c.empty()) {

        //Give up if previous chunk done mapping
        if (mappers_[ch].chunk_mapped() && !mappers_[ch].finished()) {
            mappers_[ch].request_reset();
        }
        return false;
    }

    //Start new read if mapper inactive
    if (mappers_[ch].get_state() == Mapper::State::INACTIVE) {
        mappers_[ch].new_read(c);
        active_queue_.push_back(ch);
        return true;

    } else if (mappers_[ch].get_read().number_ == c.get_number()) {

        //Don't add if previous chunk is still mapping
        if (!mappers_[ch].chunk_mapped()) {
            return false;
        }

        return mappers_[ch].add_chunk(c);
    }

    return false;
}

//TODO: make sure update is the same
std::vector<MapResult> RealtimePool::update() {

    std::vector< u16 > read_counts(threads_.size(), 0);
    u16 active_count = 0;
    std::vector<MapResult> ret;

    //Get alignment outputs
    for (u16 t = 0; t < threads_.size(); t++) {
        if (!threads_[t].out_chs_.empty()) {

            //Store and empty thread output buffer
            threads_[t].out_mtx_.lock();
            out_chs_.swap(threads_[t].out_chs_);
            threads_[t].out_mtx_.unlock();

            //Loop over alignments
            for (auto ch : out_chs_) {
                ReadBuffer &r = mappers_[ch].get_read();
                ret.emplace_back(r.get_channel(), r.number_, r.loc_);

                //TODO rename set_inactive?
                mappers_[ch].deactivate();
            }
            out_chs_.clear();
        }

        //Count reads aligning in each thread
        read_counts[t] = threads_[t].read_count();
        active_count += read_counts[t];
    }

    //Buffer queue should be ordered in "ord" mode
    for (u16 i = buffer_queue_.size()-1; i < buffer_queue_.size(); i--) {
        u16 ch = buffer_queue_[i];//TODO: store chunks in queue
        Chunk &c = chunk_buffer_[ch];

        bool added = false;

        //std::cout << "# BUFFER?\n";

        if (mappers_[ch].get_state() == Mapper::State::INACTIVE) {
            mappers_[ch].new_read(c);
            active_queue_.push_back(ch);
            added = true;
        } else if (!mappers_[ch].finished()) {
            added = mappers_[ch].add_chunk(c);
        }

        if (added) {
            if (i != buffer_queue_.size()-1) {
                buffer_queue_[i] = buffer_queue_.back();
            }
            buffer_queue_.pop_back();
        }
    }

    if (time_.get() >= 1000 && active_count > 0) {
        std::cout << "#prefill_threads ("
                  << active_count << ")";
        for (u16 c : read_counts) std::cout << " " << c;
        std::cout << "\n";
        std::cout.flush();
    }

    //Estimate how much to fill each thread
    u16 target = min(active_queue_.size() + active_count, PRMS.max_active_reads),
        min_per_thread = target / threads_.size(), // + (target % threads_.size() > 0);
        remain = target % threads_.size();

    for (u32 c : read_counts) remain -= (c > min_per_thread);
    

    u16 st = (u16) rand();

    for (u16 i = 0; i < threads_.size(); i++) {
        if (active_queue_.empty()) break;

        u16 t = (st+i) % threads_.size();

        u32 n = min_per_thread + (remain > 0) - read_counts[t];

        //If thread not full
        if (n > 0 && n <= active_queue_.size()) {
            u32 r0 = active_queue_.size() - n,
                rn = active_queue_.size();

            threads_[t].in_mtx_.lock();

            std::vector<u16> &in_chs = threads_[t].in_chs_;
            in_chs.insert(in_chs.end(),
                          &(active_queue_[r0]),
                          &(active_queue_[rn]));

            threads_[t].in_mtx_.unlock();

            active_queue_.resize(r0);
            remain -= (remain > 0);

            read_counts[t] += n;
            active_count += n;
        }
    }

    if (time_.get() >= 1000 && active_count > 0) {
        time_.reset();

        std::cout << "#pstfill_threads ("
                  << active_count << ")";

        for (u16 c : read_counts) std::cout << " " << c;
        std::cout << "\n";
        std::cout.flush();
    }

    return ret;
}

//void u32 ReadBuffer::end_read(u16 ch, u32 number) {
//    ch--;
//    if (!mappers_[ch].finished() && mappers_[ch].get_read()
//}

bool RealtimePool::all_finished() {
    if (!buffer_queue_.empty()) return false;

    for (MapperThread &t : threads_) {
        if (t.read_count() > 0 || !t.out_chs_.empty()) return false;
    }

    return true;
}

void RealtimePool::stop_all() {
    for (MapperThread &t : threads_) {
        t.running_ = false;
        t.thread_.join();
    }
}

u16 RealtimePool::MapperThread::num_threads = 0;

RealtimePool::MapperThread::MapperThread(std::vector<Mapper> &mappers)
    : tid_(num_threads++),
      mappers_(mappers),
      running_(true) {}

RealtimePool::MapperThread::MapperThread(MapperThread &&mt) 
    : tid_(mt.tid_),
      mappers_(mt.mappers_),
      running_(mt.running_),                                             
      thread_(std::move(mt.thread_)) {}

void RealtimePool::MapperThread::start() {
    thread_ = std::thread(&RealtimePool::MapperThread::run, this);
}


u16 RealtimePool::MapperThread::read_count() const {
    return in_chs_.size() + active_chs_.size();
}

void RealtimePool::MapperThread::run() {
    std::string fast5_id;
    std::vector<float> fast5_signal;

    std::vector<u16> finished;

    while (running_) {
        if (read_count() == 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            continue;
        }


        //Read inputs (pop, lock, and swap it)
        if (!in_chs_.empty()) {
            in_mtx_.lock();
            in_tmp_.swap(in_chs_);
            in_mtx_.unlock();

            for (auto ch : in_tmp_) {
                active_chs_.push_back(ch);
            }

            in_tmp_.clear(); //(pop)
        }

        //TODO: reads are in here
        //Map chunks
        for (u16 i = 0; i < active_chs_.size() && running_; i++) {
            u16 ch = active_chs_[i];

            //debug_out_ << "# thread "
            //           << tid_ << " has ch"
            //           << (ch+1) << "\n";

            mappers_[ch].process_chunk();

            mappers_[ch].thread_accs_.push_back(tid_);
            if (mappers_[ch].map_chunk()) {
                out_tmp_.push_back(i);
            }
        }

        //Add finished to output
        if (!out_tmp_.empty()) {
            out_mtx_.lock();
            for (auto i : out_tmp_) {
                out_chs_.push_back(active_chs_[i]);
            }
            out_mtx_.unlock();

            std::sort(out_tmp_.begin(), out_tmp_.end(),
                      [](u32 a, u32 b) { return a > b; });
            for (auto i : out_tmp_) {
                active_chs_[i] = active_chs_.back();
                active_chs_.pop_back();
            }
            out_tmp_.clear();
        }
    }
}
