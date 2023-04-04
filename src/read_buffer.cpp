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
#include "intervals.hpp"

//ReadBuffer::Params ReadBuffer::PRMS = {
//    flowcell : "",
//    kit : "",
//    num_channels : 512,
//    bp_per_sec   : 450,
//    sample_rate  : 4000,
//    seq_fwd      : true,
//};

ReadBuffer::ReadBuffer() {
}


ReadBuffer::ReadBuffer(const std::string &_id, u16 _channel, u32 _number, u64 start_time, 
      const std::vector<float> &raw_data, u32 raw_st, u32 raw_len) 
    : id(_id),
      channel(_channel),
      number(_number),
      start_sample(start_time) {
    //signal = ValArray<float>(&raw_data[raw_st], raw_len);
    set_signal(&(raw_data[raw_st]), raw_len);
    //if (raw_st + raw_len > raw_data.size()) raw_len = raw_data.size() - raw_st;
    //signal.resize(raw_len);
    //for (u32 i = 0; i < raw_len; i++) signal[i] = raw_data[raw_st+i];
}

void ReadBuffer::set_signal(const float *ptr, size_t len) {
    signal = ValArray<float>(ptr, len);
}

#ifdef PYBIND
ReadBuffer::ReadBuffer(const std::string &_id, u16 _channel, u32 _number, u64 start_time, const py::array_t<float> &raw_data) 
    : id(_id),
      channel(_channel),
      number(_number),
      start_sample(start_time) {

    PyArray<float> arr(raw_data);
    set_signal(arr.data, arr.size());
    //signal.assign(arr.begin(), arr.end());
    //signal.reserve(arr.size());
    //for (auto s : arr) {
    //    signal.push_back(s);
    //}
}
#endif

u64 ReadBuffer::size() const {
    return signal.size();
}

u64 ReadBuffer::get_start() const {
    return start_sample;
}  

u64 ReadBuffer::get_end() const {
    return start_sample + size();
}

//Returns sample at specified index
float &ReadBuffer::operator[](u32 i) {
    return signal[i];
}

bool ReadBuffer::empty() const {
    return signal.empty();
}

void ReadBuffer::clear() {
    set_signal(NULL, 0);
}

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2) {
    return r1.start_sample < r2.start_sample;
}

