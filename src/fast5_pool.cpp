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
#include "fast5_reader.hpp"
#include "fast5_pool.hpp"
#include "mapper.hpp"

bool open_fast5(const std::string &filename, fast5::File &file) {
    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" << filename << "'\n";
            return false;
        }

        return true;
        
    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        return false;
    }


    return false;
}


Fast5Pool::Fast5Pool(MapperParams &params, u16 nthreads, u32 batch_size) {
    nthreads_ = nthreads;
    batch_size_ = batch_size;

    for (u16 i = 0; i < nthreads_; i++) {
        //threads_.push_back(MapperThread(params));
        threads_.emplace_back(params);
    }
    for (MapperThread &t : threads_) {
        t.start();
    }
}

void Fast5Pool::add_fast5s(const std::vector<std::string> &new_fast5s) {
    fast5s_.insert(fast5s_.end(), new_fast5s.begin(), new_fast5s.end());
}

std::vector<std::string> Fast5Pool::update() {
    std::vector<std::string> ret;

    for (MapperThread &t : threads_) {
        t.out_mtx_.lock();
        while (!t.locs_out_.empty()) {
            //TODO: parse outside
            ret.push_back(t.locs_out_.front().str());
            t.locs_out_.pop_front();
        }
        t.out_mtx_.unlock();

        if (t.signals_in_.size() < batch_size_ / 2) {
            while (!fast5s_.empty() && t.signals_in_.size() < batch_size_) {

                std::string fname = fast5s_.front();
                fast5s_.pop_front();

                //fast5::File fast5;
                //open_fast5(fname, fast5);            
                //std::vector<float> samples = fast5.get_raw_samples();
                //std::string ids = fast5.get_raw_samples_params().read_id;
                //t.in_mtx_.lock();
                //t.signals_in_.push_back(samples);
                //t.ids_in_.push_back(ids);
                //fast5.close();
                //t.in_mtx_.unlock();

                std::vector<Fast5Read> reads;
                load_multi_fast5(fname, reads);
                t.in_mtx_.lock();
                for (Fast5Read &r : reads) {
                    t.signals_in_.push_back(r.raw_data);
                    t.ids_in_.push_back(r.id);
                }
                t.in_mtx_.unlock();
            }
        }
    }
    
    return ret;
}

bool Fast5Pool::all_finished() {
    if (!fast5s_.empty()) return false;

    for (MapperThread &t : threads_) {
        if (t.aligning_ || !t.locs_out_.empty()) return false;
    }

    return true;
}

void Fast5Pool::stop_all() {
    for (MapperThread &t : threads_) {
        t.running_ = false;
        t.thread_.join();
    }
}


Fast5Pool::MapperThread::MapperThread(MapperParams &params)
    : running_(true),
      aligning_(false),
      mapper_(params, 0) {}

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
