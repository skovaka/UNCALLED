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
#include "channel_pool.hpp"
#include "mapper.hpp"
#include "fast5_reader.hpp"

ChannelPool::ChannelPool(MapperParams &params, u16 nthreads, u16 nchannels) {
    fast5_count_ = 0;

    for (u16 t = 0; t < nthreads; t++) {
        threads_.emplace_back(mappers_);
        thread_ids_.push_back(t);
    }
    
    //mappers_.reserve(nchannels);
    channel_busy_.reserve(nchannels);
    read_buffer_.resize(nchannels);
    channel_fast5s_.resize(nchannels);
    for (u16 i = 0; i < nchannels; i++) {
        mappers_.push_back(Mapper(params, i));
        channel_busy_.push_back(false);
    }

    for (u16 t = 0; t < nthreads; t++) {
        threads_[t].start();
    }
}

void ChannelPool::add_fast5s(const std::vector<std::string> &fast5s, 
                             const std::vector<u16> &channels) {

    if (fast5s.size() != channels.size()) {
        std::cerr << "Error: size of filenames and channels not same\n";
        return;
    }

    for (u64 i = 0; i < fast5s.size(); i++) {
        channel_fast5s_[channels[i]-1].push_back(fast5s[i]);
        fast5_count_++;
    }

    for (u16 i = 0; i < channel_fast5s_.size(); i++) {
        if (channel_fast5s_[i].empty()) continue;
        read_buffer_[i] = Fast5Read(channel_fast5s_[i].front());
        channel_fast5s_[i].pop_front();
        channel_queue_.push_back(i);
        fast5_count_--;
    }

    std::cerr << "Reads " << fast5_count_ << " " << channel_queue_.size() << "\n";
}

std::vector<std::string> ChannelPool::update() {
    std::vector<std::string> ret;

    std::vector< u16 > read_counts(threads_.size());
    u16 active_count = 0;

    //Get alignment outputs
    std::deque<ReadLoc> locs;
    for (u16 t = 0; t < threads_.size(); t++) {
        if (!threads_[t].outputs_.empty()) {

            //Store and empty thread output buffer
            threads_[t].out_mtx_.lock();
            locs.swap(threads_[t].outputs_);
            threads_[t].out_mtx_.unlock();

            //Loop over alignments
            while(!locs.empty()) {

                //Mark channel as not aligning
                u16 ch = locs.front().get_channel();
                channel_busy_[ch] = false;

                //If channel has next read, add channel to aln queue
                if (!read_buffer_[ch].empty()) {
                    channel_queue_.push_back(ch);
                }

                //Store previous alignment
                ret.push_back(locs.front().str());
                locs.pop_front();
            }
        }

        //Count reads aligning in each thread
        read_counts[t] = threads_[t].read_count();
        active_count += read_counts[t];
    }

    //Estimate how much to fill each thread
    u16 target = channel_queue_.size() + active_count,
        per_thread = target / threads_.size() + (target % thread_ids_.size() > 0);

    for (u16 t = 0; t < thread_ids_.size(); t++) {

        //If thread not full
        if (read_counts[t] < per_thread) {
            threads_[t].in_mtx_.lock();

            //Fill thread till full
            while (!channel_queue_.empty() && read_counts[t] < per_thread) {

                //Get next channel to align
                u16 ch = channel_queue_.front(); 
                channel_queue_.pop_front();

                //Add new read to thread queue
                threads_[t].new_reads_.push_back(Fast5Read());
                threads_[t].new_reads_.back().swap(read_buffer_[ch]);

                //Mark channel busy
                channel_busy_[ch] = true;
                read_counts[t]++;
            }
            threads_[t].in_mtx_.unlock();
        }
    }

    Timer t;

    //Open new fast5s
    for (u16 ch = 0; ch < channel_fast5s_.size(); ch++) {
        if (!read_buffer_[ch].empty() || channel_fast5s_[ch].empty()) continue;
        read_buffer_[ch] = Fast5Read(channel_fast5s_[ch].front());
        channel_fast5s_[ch].pop_front();
        fast5_count_--;
        if (!channel_busy_[ch]) channel_queue_.push_back(ch);
    }

    return ret;
}

bool ChannelPool::all_finished() {
    if (fast5_count_ > 0) return false;

    for (MapperThread &t : threads_) {
        if (t.read_count() > 0 || !t.outputs_.empty()) return false;
    }

    return true;
}

void ChannelPool::stop_all() {
    for (MapperThread &t : threads_) {
        t.running_ = false;
        t.thread_.join();
    }
}

u16 ChannelPool::MapperThread::num_threads = 0;

ChannelPool::MapperThread::MapperThread(std::vector<Mapper> &mappers)
    : tid_(num_threads++),
      mappers_(mappers),
      running_(true) {}

ChannelPool::MapperThread::MapperThread(MapperThread &&mt) 
    : tid_(mt.tid_),
      mappers_(mt.mappers_),
      running_(mt.running_),                                             
      thread_(std::move(mt.thread_)) {}

void ChannelPool::MapperThread::start() {
    thread_ = std::thread(&ChannelPool::MapperThread::run, this);
}

u16 ChannelPool::MapperThread::read_count() const {
    return new_reads_.size() + active_reads_.size();
}

void ChannelPool::MapperThread::run() {
    std::string fast5_id;
    std::vector<float> fast5_signal;

    std::vector<u16> finished;

    while (running_) {
        if (read_count() == 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }

        if (!new_reads_.empty()) {
            //std::cerr << "New read " << tid_ << "\n";
            std::vector<Fast5Read> new_tmp;
            in_mtx_.lock();
            new_reads_.swap(new_tmp);
            in_mtx_.unlock();

            for (Fast5Read &r : new_tmp) {
                mappers_[r.channel].new_read(r.id, r.number);
                active_reads_.push_back(Fast5Read());
                active_reads_.back().swap(r);
            }
        }


        for (u16 i = 0; i < active_reads_.size() && running_; i++) {
            //std::cerr << "Aligning " << tid_ << " " << i << "\n";
            std::vector<Fast5Read> new_tmp;
            Fast5Read &r = active_reads_[i];
            if (r.empty() || mappers_[r.channel].add_sample(r.next_sig())) {
                finished.push_back(i);
            }
        }

        if (!finished.empty()) {
            out_mtx_.lock();
            for (u16 f : finished) {
                outputs_.push_back(mappers_[active_reads_[f].channel].get_loc());
            }
            out_mtx_.unlock();

            while(!finished.empty()) {
                u16 f = finished.back();
                finished.pop_back();
                active_reads_[f].swap(active_reads_.back());
                active_reads_.pop_back();
            }
        }
    }
}
