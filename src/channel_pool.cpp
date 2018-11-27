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


Fast5Read::Fast5Read() :
    signal(),
    name(""),
    channel(0),
    i(0) {}

Fast5Read::Fast5Read(const std::string &filename) 
    : i(0) {

    if (filename.empty()) {
        signal = std::vector<float>();
        name = "";
        channel = 0;
        return;
    }

    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
        return;
    }

    fast5::File file;

    try {
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" << filename << "'\n";
        } else {
            
            signal = file.get_raw_samples();
            name = file.get_raw_samples_params().read_id;
            channel = atoi(file.get_channel_id_params().channel_number.c_str()) - 1;
        
        }

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }
}

void Fast5Read::swap(Fast5Read &r) {
    signal.swap(r.signal);
    name.swap(r.name);
    std::swap(i, r.i);
    std::swap(channel, r.channel);
}

float Fast5Read::next_sig() {
    return signal[i++];
}

bool Fast5Read::empty() const {
    return i >= signal.size();
}
    

ChannelPool::ChannelPool(MapperParams &params, u16 nthreads, u16 nchannels) {
    
    for (u16 t = 0; t < nthreads; t++) {
        threads_.emplace_back(mappers_);
        thread_ids_.push_back(t);
    }
    
    //mappers_.reserve(nchannels);
    channel_busy_.reserve(nchannels);
    read_buffer_.resize(nchannels);
    for (u16 i = 0; i < nchannels; i++) {
        mappers_.push_back(Mapper(params, i));
        channel_busy_.push_back(false);
    }

    for (u16 t = 0; t < nthreads; t++) {
        threads_[t].start();
    }
}

void ChannelPool::add_fast5s(const std::vector<std::string> &new_fast5s) {
    filenames_.insert(filenames_.end(), new_fast5s.begin(), new_fast5s.end());

    if (next_read_.empty()) {
        next_read_ = Fast5Read(filenames_.front());
        filenames_.pop_front();
        //std::cout << next_read_.channel << " init\n";
    }
}

std::vector<std::string> ChannelPool::update() {
    std::deque<ReadLoc> locs;
    std::vector<std::string> ret;

    //Fill channel buffers
    while (!next_read_.empty() && read_buffer_[next_read_.channel].empty()) {

        //If no other read from that channel is being aligned, queue up that channel
        //Only necissary for fast5 alignment - will not happen iin real-time
        if (!channel_busy_[next_read_.channel]) {
            channel_queue_.push_back(next_read_.channel);
        }

        //Add next read to buffer
        read_buffer_[next_read_.channel].swap(next_read_);

        //Load next read from filename list
        if (filenames_.empty()) break;

        next_read_ = Fast5Read(filenames_.front());
        filenames_.pop_front();
    }

    std::vector< u16 > in_counts(threads_.size());

    //Get outputs
    for (u16 t = 0; t < threads_.size(); t++) {
        threads_[t].out_mtx_.lock();
        locs.swap(threads_[t].outputs_);
        threads_[t].out_mtx_.unlock();

        while(!locs.empty()) {
            u16 ch = locs.front().get_channel();
            channel_busy_[ch] = false;
            if (!read_buffer_[ch].empty()) {
                channel_queue_.push_back(ch);
            }
            ret.push_back(locs.front().str());
            locs.pop_front();
        }

        in_counts[t] = threads_[t].inputs_.size();
    }

    std::sort(thread_ids_.begin(), thread_ids_.end(),
                [&in_counts](u16 t1, u16 t2) {
                    return in_counts[t1] < in_counts[t2];
                });

    for (u16 i = 0; i < thread_ids_.size(); i++) {
        MapperThread &t = threads_[thread_ids_[i]];
        u16 next_count = in_counts[thread_ids_[(i+1)%threads_.size()]];

        t.in_mtx_.lock();
        while(t.inputs_.size() <= next_count && !channel_queue_.empty()) {
            u16 ch = channel_queue_.front();
            t.inputs_.push_back(Fast5Read());
            t.inputs_.back().swap(read_buffer_[ch]);
            channel_busy_[ch] = true;
            channel_queue_.pop_front();
        }
        t.in_mtx_.unlock();
    }
        
    return ret;
}

bool ChannelPool::all_finished() {
    if (!filenames_.empty()) return false;

    for (MapperThread &t : threads_) {
        if (!t.inputs_.empty() || !t.outputs_.empty()) return false;
    }

    return true;
}

void ChannelPool::stop_all() {
    for (MapperThread &t : threads_) {
        t.running_ = false;
        t.thread_.join();
    }
}


ChannelPool::MapperThread::MapperThread(std::vector<Mapper> &mappers)
    : mappers_(mappers),
      running_(true) {}

ChannelPool::MapperThread::MapperThread(MapperThread &&mt) 
    : mappers_(mt.mappers_),
      running_(mt.running_),                                             
      thread_(std::move(mt.thread_)) {}

void ChannelPool::MapperThread::start() {
    thread_ = std::thread(&ChannelPool::MapperThread::run, this);
}

void ChannelPool::MapperThread::run() {
    std::string fast5_id;
    std::vector<float> fast5_signal;

    std::vector<u16> finished;

    while (running_) {
        if (inputs_.empty()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }

        in_mtx_.lock();
        u16 num_reads = inputs_.size();
        in_mtx_.unlock();

        for (u16 i = 0; i < num_reads; i++) {

            Fast5Read &r = inputs_[i];
            //TODO: don't check every time?
            if (r.i == 0) {
                mappers_[r.channel].new_read(r.name);
            }

            //TODO: add signal in chunks?
            if (r.empty() || mappers_[r.channel].add_sample(r.next_sig())) {
                finished.push_back(i);
            }
        }

        if (!finished.empty()) {
            out_mtx_.lock();
            for (u16 f : finished) {
                outputs_.push_back(mappers_[inputs_[f].channel].get_mapping());
            }
            out_mtx_.unlock();

            in_mtx_.lock();
            while(!finished.empty()) {
                u16 f = finished.back();
                finished.pop_back();
                inputs_[f].swap(inputs_.back());
                inputs_.pop_back();
            }
            in_mtx_.unlock();
        }

    }
}
