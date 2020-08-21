#ifndef SYNC_OUT_HPP
#define SYNC_OUT_HPP

#include <ostream>
#include <mutex>



#define SYNC_OUT(o, v) \
    sync_out_mtx_.lock(); \
    o << v; \
    sync_out_mtx_.unlock(); 

class SyncOut {
    public:
    static std::mutex mtx_;

    std::ostream &out_;

    public:
    SyncOut(std::ostream &out) : out_ (out) {}
    void write(const std::string &s);
};

#endif
