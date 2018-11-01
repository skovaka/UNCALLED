#ifndef MAP_POOL_HPP
#define MAP_POOL_HPP

#include <thread>
#include <vector>
#include <list>
#include "mapper.hpp"

class Fast5Pool {
    public:
    Fast5Pool(MapperParams &params, u16 nthreads, u32 batch_size);
    
    void add_fast5s(const std::list<std::string> &new_fast5s);
    std::vector<std::string> update();
    bool all_finished();
    void stop_all();

    private:

    class MapperThread {
        public:
        MapperThread(MapperParams &params);
        MapperThread(MapperThread &&mt);

        void start();
        void run();

        bool running_, aligning_;
        Mapper mapper_;
        std::thread thread_;
        std::list<std::string> ids_in_;
        std::list< std::vector<float> > signals_in_;
        std::list<ReadLoc> locs_out_;
        std::mutex in_mtx_, out_mtx_;
        //std::unique_lock<std::mutex> in_lck_, out_lck_;
    };

    u16 nthreads_;
    std::vector<MapperThread> threads_;
    std::list<std::string> fast5s_;
    u32 batch_size_;
};


#endif
