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

#include <iostream>
#include "read_buffer.hpp"

ReadBuffer::Params ReadBuffer::PRMS = {
    num_channels : 512,
    bp_per_sec   : 450,
    sample_rate  : 4000,
    chunk_time   : 1.0,
    start_chunk  : 0,
    max_chunks   : 1000000,
    seq_fwd      : true,
    skip_notempl : false
};

ReadBuffer::ReadBuffer() {
}

ReadBuffer::ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, const std::string &dtype, const std::string &raw_str)
    : id_(id),
      channel_idx_(channel-1),
      number_(number),
      start_sample_(start_time),
      chunk_count_(1) {

    //TODO: could store chunk data as C arrays to prevent extra copy
    //probably not worth it
    if (dtype == "float32") {
        signal_.resize(raw_str.size()/sizeof(float));
        float *raw_arr = (float *) raw_str.data();
        signal_.assign(raw_arr, &raw_arr[signal_.size()]);

    } else if (dtype == "int16") {
        signal_.resize(raw_str.size()/sizeof(u16));
        i16 *raw_arr = (i16 *) raw_str.data();
        signal_.assign(raw_arr, &raw_arr[signal_.size()]);
        //for (u32 i = 0; i < signal_.size(); i++) {
        //    signal_[i] = ReadBuffer::calibrate(get_channel(), raw_arr[i]);
        //}

    } else if (dtype == "int32") {
        signal_.resize(raw_str.size()/sizeof(u32));
        i32 *raw_arr = (i32 *) raw_str.data();
        signal_.assign(raw_arr, &raw_arr[signal_.size()]);
        //for (u32 i = 0; i < signal_.size(); i++) {
        //    signal_[i] = ReadBuffer::calibrate(get_channel(), raw_arr[i]);
        //}

    } else {
        std::cerr << "Error: unsuportted raw signal dtype\n";
    }

    full_duration_ = signal_.size();
}

ReadBuffer::ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, 
      const std::vector<float> &raw_data, u32 raw_st, u32 raw_len) 
    : id_(id),
      channel_idx_(channel-1),
      number_(number),
      start_sample_(start_time),
      full_duration_(raw_len),
      chunk_count_(1) {
    if (raw_st + raw_len > raw_data.size()) raw_len = raw_data.size() - raw_st;
    signal_.resize(raw_len);
    for (u32 i = 0; i < raw_len; i++) signal_[i] = raw_data[raw_st+i];
}


//Returns read ID (name)
std::string ReadBuffer::get_id() const {
    return id_;
}

u16 ReadBuffer::get_channel() const {
    return channel_idx_+1;
}

u16 ReadBuffer::get_channel_idx() const {
    return channel_idx_;
}

void ReadBuffer::set_channel(u16 ch) {
    channel_idx_ = ch-1;
}

u32 ReadBuffer::get_number() const {
    return number_;
}

u64 ReadBuffer::size() const {
    return signal_.size();
}

u64 ReadBuffer::get_full_duration() const {
    return full_duration_;
}

void ReadBuffer::set_full_duration(u64 raw_len) {
    full_duration_ = raw_len;
}

u64 ReadBuffer::get_start() const {
    return start_sample_;
}  

void ReadBuffer::set_start(u64 start) {
    start_sample_ = start;
}  

u64 ReadBuffer::get_end() const {
    return start_sample_ + full_duration_;
}

//Returns the signal buffer
const std::vector<float> &ReadBuffer::get_signal() const {
    return signal_;
}

//Returns sample at specified index
float &ReadBuffer::operator[](u32 i) {
    return signal_[i];
}

bool ReadBuffer::empty() const {
    return signal_.empty();
}

void ReadBuffer::clear() {
    signal_.clear();
}

bool ReadBuffer::chunks_maxed() const {
    return chunk_count_ >= PRMS.max_chunks;
}

u32 ReadBuffer::get_chunk_count() const {
    return chunk_count_; //signal_.size() / PRMS.chunk_len();
}

ReadBuffer ReadBuffer::get_chunk(u32 i) const {
    u32 st = i * PRMS.chunk_len(),
        ln = PRMS.chunk_len();

    if (st > signal_.size()) {
        st = signal_.size();
    } 
    
    if (st+ln > signal_.size()) {
        ln = signal_.size() - st;
    }

    return ReadBuffer(id_, get_channel(), number_, start_sample_+st, signal_, st, ln);
}

//TODO eliminate
u32 ReadBuffer::get_chunks(std::vector<ReadBuffer> &chunks, bool real_start, u32 offs) const {
    u32 count = 0;
    u16 l = PRMS.chunk_len();

    float start = real_start ? start_sample_ : 0;

    for (u32 i = offs; i+l <= signal_.size() && count < PRMS.max_chunks; i += l) {
        chunks.emplace_back(id_, get_channel(), number_, 
                                 start+i, signal_, i, l);
        count++;
    }
    return count;

}

bool ReadBuffer::add_next_chunk(ReadBuffer &c) {
    if (channel_idx_ != c.get_channel_idx() || 
        number_ != c.get_number()) return false;

    signal_.swap(c.signal_);
    c.signal_.clear();
    full_duration_ += signal_.size();
    chunk_count_++;

    return true;
}

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2) {
    return r1.start_sample_ < r2.start_sample_;
}

