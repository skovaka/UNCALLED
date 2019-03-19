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
#include "chunk_pool.hpp"
#include "mapper.hpp"

//TODO: legit chunk default
ReadStream::ReadStream(u16 chunk_len=4500, const std::string &name="") {
    name_ = name;
}

void ReadStream::set_name(const std::string &name) {
    name_ = name;
}

bool ReadStream::empty() const {
    return chunk_.empty();
}

ChunkPool::ChunkPool(MapperParams &params, u16 nthreads, u16 nchannels) {
    for (u16 t = 0; t < nthreads; t++) {
        threads_.emplace_back(mappers_);
    }
    
    //mappers_.reserve(nchannels);
    channel_active_.reserve(nchannels);
    read_buffer_.resize(nchannels);
    buffer_queue_.reserve(nchannels);
    for (u16 i = 0; i < nchannels; i++) {
        mappers_.push_back(Mapper(params, i));
        channel_active_.push_back(false);
    }

    for (u16 t = 0; t < nthreads; t++) {
        threads_[t].start();
    }
}

void ChunkPool::new_read(u16 ch, const std::string &name) {
    read_buffer_[ch].name_ = name;
    //TODO: what if channel currently active?
    //in practice rare, but could happen, especially in naive simulation
}

//Add chunk to master buffer
bool ChunkPool::add_chunk(u16 ch, std::vector<float> &chunk) {
    if (!read_buffer_[ch].empty()) return false; //TODO: reset mapper here?
    chunk.swap(read_buffer_[ch].chunk_);
    buffer_queue_.push_back(ch);
    return true;
}

std::vector<std::string> ChunkPool::update() {
    std::vector<std::string> ret;

    std::vector< u16 > read_counts(threads_.size());
    u16 active_count = 0;

    //Get alignment outputs
    //TODO: redo this
    for (u16 t = 0; t < threads_.size(); t++) {
        if (!threads_[t].out_chs_.empty()) {

            //Store and empty thread output buffer
            threads_[t].out_mtx_.lock();
            out_chs_.swap(threads_[t].out_chs_);
            threads_[t].out_mtx_.unlock();

            //Loop over alignments
            for (auto ch : out_chs_) {
                channel_active_[ch] = false;
                ret.push_back(mappers_[ch].get_loc().str());
            }
            out_chs_.clear();
        }

        //Count reads aligning in each thread
        read_counts[t] = threads_[t].read_count();
        active_count += read_counts[t];
    }
    
    //TODO: combine queues, partition into active/buffer regions
    for (u16 i = buffer_queue_.size()-1; i < buffer_queue_.size(); i--) {
        u16 ch = buffer_queue_[i];
        if (mappers_[ch].swap_chunk(read_buffer_[ch].chunk_)) {
            if (i != buffer_queue_.size()-1) {
                 buffer_queue_[i] = buffer_queue_.back();
            }
            buffer_queue_.pop_back();
            if (!channel_active_[ch]) {
                active_queue_.push_back(ch);
                mappers_[ch].new_read(read_buffer_[ch].name_);
            }
        } 
    }

    //Estimate how much to fill each thread
    u16 target = active_queue_.size() + active_count,
        per_thread = target / threads_.size() + (target % threads_.size() > 0);

    for (u16 t = 0; t < threads_.size(); t++) {

        //If thread not full
        if (read_counts[t] < per_thread) {

            //Fill thread till full
            //TODO: compute number exactly, only lock while adding
            threads_[t].in_mtx_.lock();
            while (!active_queue_.empty() && read_counts[t] < per_thread) {
                u16 ch = active_queue_.back(); 
                active_queue_.pop_back();
                threads_[t].in_chs_.push_back(ch);
                channel_active_[ch] = true;
                read_counts[t]++;
            }
            threads_[t].in_mtx_.unlock();
        }
    }

    return ret;
}

bool ChunkPool::all_finished() {
    if (!buffer_queue_.empty()) return false;

    for (MapperThread &t : threads_) {
        if (t.read_count() > 0 || !t.out_chs_.empty()) return false;
    }

    return true;
}

void ChunkPool::stop_all() {
    for (MapperThread &t : threads_) {
        t.running_ = false;
        t.thread_.join();
    }
}

u16 ChunkPool::MapperThread::num_threads = 0;

ChunkPool::MapperThread::MapperThread(std::vector<Mapper> &mappers)
    : tid_(num_threads++),
      mappers_(mappers),
      running_(true) {}

ChunkPool::MapperThread::MapperThread(MapperThread &&mt) 
    : tid_(mt.tid_),
      mappers_(mt.mappers_),
      running_(mt.running_),                                             
      thread_(std::move(mt.thread_)) {}

void ChunkPool::MapperThread::start() {
    thread_ = std::thread(&ChunkPool::MapperThread::run, this);
}


u16 ChunkPool::MapperThread::read_count() const {
    return in_chs_.size() + active_chs_.size();
}

void ChunkPool::MapperThread::run() {
    std::string fast5_id;
    std::vector<float> fast5_signal;

    std::vector<u16> finished;

    while (running_) {
        if (read_count() == 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }

        //Read inputs (pop, lock, and swap it)
        if (!in_chs_.empty()) {
            in_mtx_.lock();
            in_tmp_.swap(in_chs_);
            in_mtx_.unlock();

            for (auto ch : in_tmp_) {
                mappers_[ch].process_chunk();
                active_chs_.push_back(ch);
            }

            in_tmp_.clear(); //(pop)
        }

        //Map chunks
        for (u16 i = 0; i < active_chs_.size() && running_; i++) {
            u16 ch = active_chs_[i];
            if (mappers_[ch].map_chunk()) {
                out_tmp_.push_back(i);
            }
        }

        //Add finished to output
        if (!out_tmp_.empty()) {
            out_mtx_.lock();
            for (auto i : out_tmp_) out_chs_.push_back(active_chs_[i]);
            out_mtx_.unlock();

            for (auto i : out_tmp_) {
                active_chs_[i] = active_chs_.back();
                active_chs_.pop_back();
            }
            out_tmp_.clear();
        }
    }
}
