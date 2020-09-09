#include <iostream>
#include "chunk.hpp"
#include "read_buffer.hpp"

std::vector<float> Chunk::cal_offsets_,
                   Chunk::cal_coefs_;

Chunk::Chunk() 
    : id_(""),
      channel_idx_(0),
      number_(0),
      start_time_(0),
      raw_data_() {}


Chunk::Chunk(const std::string &id, u16 channel, u32 number, u64 chunk_start, 
             const std::string &dtype, const std::string &raw_str) 
    : id_(id),
      channel_idx_(channel-1),
      number_(number),
      start_time_(chunk_start) {

    //TODO: could store chunk data as C arrays to prevent extra copy
    //probably not worth it
    if (dtype == "float32") {
        raw_data_.resize(raw_str.size()/sizeof(float));
        float *raw_arr = (float *) raw_str.data();
        raw_data_.assign(raw_arr, &raw_arr[raw_data_.size()]);

    } else if (dtype == "int16") {
        raw_data_.resize(raw_str.size()/sizeof(u16));
        i16 *raw_arr = (i16 *) raw_str.data();
        raw_data_.assign(raw_arr, &raw_arr[raw_data_.size()]);
        //for (u32 i = 0; i < raw_data_.size(); i++) {
        //    raw_data_[i] = ReadBuffer::calibrate(get_channel(), raw_arr[i]);
        //}

    } else if (dtype == "int32") {
        raw_data_.resize(raw_str.size()/sizeof(u32));
        i32 *raw_arr = (i32 *) raw_str.data();
        raw_data_.assign(raw_arr, &raw_arr[raw_data_.size()]);
        //for (u32 i = 0; i < raw_data_.size(); i++) {
        //    raw_data_[i] = ReadBuffer::calibrate(get_channel(), raw_arr[i]);
        //}

    } else {
        std::cerr << "Error: unsuportted raw signal dtype\n";
    }

    //TODO: templatize
}

Chunk::Chunk(const std::string &id, u16 channel, u32 number, u64 start_time, 
             const std::vector<float> &raw_data, u32 raw_st, u32 raw_len) 
    : id_(id),
      channel_idx_(channel-1),
      number_(number),
      start_time_(start_time) {
    if (raw_st + raw_len > raw_data.size()) raw_len = raw_data.size() - raw_st;
    raw_data_.resize(raw_len);
    for (u32 i = 0; i < raw_len; i++) raw_data_[i] = raw_data[raw_st+i];
    
}
//Chunk::Chunk(const Chunk &c) 
//    : id_(c.id_),
//      channel_idx_(c.channel_idx_),
//      number_(c.number_),
//      start_time_(c.start_time_),
//      raw_data_(c.raw_data_) {}

float &Chunk::operator[] (u32 i) {
    return raw_data_[i];
}

u32 Chunk::size() const {
    return raw_data_.size();
}

bool Chunk::empty() const {
    return raw_data_.empty();
}

void Chunk::print() const {
    for (float s : raw_data_) std::cout << s << std::endl;
}

bool Chunk::pop(std::vector<float> &raw_data) {
    raw_data_.swap(raw_data);
    clear();
    return !raw_data.empty();
}

u64 Chunk::get_start() const {
    return start_time_;
}

u64 Chunk::get_end() const {
    return start_time_ + raw_data_.size();
}

std::string Chunk::get_id() const {
    return id_;
}

u16 Chunk::get_channel_idx() const {
    return channel_idx_;
}

u16 Chunk::get_channel() const {
    return channel_idx_+1;
}

u32 Chunk::get_number() const {
    return number_;
}

void Chunk::set_start(u64 time) {
    start_time_ = time;
}

void Chunk::swap(Chunk &c) {
    std::swap(id_, c.id_);
    std::swap(channel_idx_, c.channel_idx_);
    std::swap(number_, c.number_);
    std::swap(start_time_, c.start_time_);
    raw_data_.swap(c.raw_data_);
}

void Chunk::clear() {
    raw_data_.clear();
}   

bool operator< (const Chunk &r1, const Chunk &r2) {
    return r1.start_time_ < r2.start_time_;
}
