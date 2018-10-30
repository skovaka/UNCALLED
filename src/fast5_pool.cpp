#include <thread>
#include "fast5_pool.hpp"
#include "mapper.hpp"


Fast5Pool::Fast5Pool(MapperParams &params, u16 nthreads, u32 batch_size) {
    nthreads_ = nthreads;
    batch_size_ = batch_size;

    for (u16 i = 0; i < nthreads_; i++) {
        threads_.push_back(MapperThread(params));
    }
    for (MapperThread &t : threads_) {
        t.start();
    }
}

void Fast5Pool::add_fast5s(const std::list<std::string> &new_fast5s) {
    fast5s_.insert(fast5s_.end(), new_fast5s.begin(), new_fast5s.end());
}

std::vector<std::string> Fast5Pool::update() {
    std::vector<std::string> ret;

    for (MapperThread &t : threads_) {
        t.out_mtx_.lock(); //maybe only nescissary for last element?
        while (!t.locs_out_.empty()) {
            ret.push_back(t.locs_out_.front().str());
            t.locs_out_.pop_front();
        }
        t.out_mtx_.unlock();

        if (t.fast5s_in_.size() < batch_size_ / 2) {
            t.in_mtx_.lock(); //maybe only nescissary for last element?
            while (!fast5s_.empty() && t.fast5s_in_.size() < batch_size_) {
                t.fast5s_in_.push_back(fast5s_.front());
                fast5s_.pop_front();
            }
            t.in_mtx_.unlock();
        }
    }
    
    return ret;
}

bool Fast5Pool::all_finished() {
    if (!fast5s_.empty()) return false;

    for (MapperThread &t : threads_) {
        if (t.aligning_) return false;
    }

    return true;
}

void Fast5Pool::stop_all() {
    for (MapperThread &t : threads_) {
        t.running_ = false;
    }
}


Fast5Pool::MapperThread::MapperThread(MapperParams &params)
    : running_(true),
      aligning_(false),
      mapper_(Mapper(params)) {}

Fast5Pool::MapperThread::MapperThread(MapperThread &&mt) 
    : running_(mt.running_),                                             
      aligning_(mt.aligning_),                                           
      mapper_(std::move(mt.mapper_)) 
{

}

void Fast5Pool::MapperThread::start() {
    thread_ = std::thread(&Fast5Pool::MapperThread::run, this);
}

void Fast5Pool::MapperThread::run() {
    std::string fast5_name;
    ReadLoc loc;
    while (running_) {
        if (fast5s_in_.empty()) {
            aligning_ = false;
            //sleep?
            continue;
        }

        aligning_ = true;
        
        fast5_name = fast5s_in_.front();
        in_mtx_.lock(); //check if last elem?
        fast5s_in_.pop_front();
        in_mtx_.unlock();

        loc = mapper_.map_fast5(fast5_name);

        out_mtx_.lock(); //you get it
        locs_out_.push_back(loc);
        out_mtx_.unlock();
    }
}
