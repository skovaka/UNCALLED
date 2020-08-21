#include "sync_out.hpp"

std::mutex SyncOut::mtx_;
SyncOut out(std::cout);

void SyncOut::write(const std::string &s) {
    mtx_.lock();
    out_.write(s);
    mtx_.unlock();
}
