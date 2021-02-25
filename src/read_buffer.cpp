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

#include "read_buffer.hpp"

ReadBuffer::Params ReadBuffer::PRMS = {
    num_channels : 512,
    bp_per_sec   : 450,
    sample_rate  : 4000,
    chunk_time   : 1.0,
    max_chunks   : 1000000,
};

ReadBuffer::ReadBuffer() {
    chunk_count_ = 0;
    
}

//TODO: eliminate swap from mapper, rely on automatic move constructor
void ReadBuffer::swap(ReadBuffer &r) {
    //std::swap(source_, r.source_);
    std::swap(channel_idx_, r.channel_idx_);
    std::swap(id_, r.id_);
    std::swap(number_, r.number_);
    std::swap(start_sample_, r.start_sample_);
    std::swap(raw_len_, r.raw_len_);
    std::swap(full_signal_, r.full_signal_);
    std::swap(chunk_, r.chunk_);
    std::swap(chunk_count_, r.chunk_count_);
    std::swap(chunk_processed_, r.chunk_processed_);
}

void ReadBuffer::clear() {
    raw_len_ = 0;
    full_signal_.clear();
    chunk_.clear();
    chunk_count_ = 0;
}

ReadBuffer::ReadBuffer(const hdf5_tools::File &file, 
                       const std::string &raw_path, 
                       const std::string &ch_path, 
                       const std::string &seg_path) {

    for (auto a : file.get_attr_map(raw_path)) {
        if (a.first == "read_id") {
            id_ = a.second;
        } else if (a.first == "read_number") {
            number_ = atoi(a.second.c_str());
        } else if (a.first == "start_time") {
            start_sample_ = atoi(a.second.c_str());
        }
    }

	float cal_digit = 1, cal_range = 1, cal_offset = 0;
    for (auto a : file.get_attr_map(ch_path)) {
        if (a.first == "channel_number") {
            channel_idx_ = atoi(a.second.c_str()) - 1;
        } else if (a.first == "digitisation") {
            cal_digit = atof(a.second.c_str());
        } else if (a.first == "range") {
            cal_range = atof(a.second.c_str());
        } else if (a.first == "offset") {
            cal_offset = atof(a.second.c_str());
        }
    }

    u32 samp_st = 0;

    if (!seg_path.empty()) {
        for (auto a : file.get_attr_map(seg_path)) {
            if (a.first == "first_sample_template") {
                samp_st = atoi(a.second.c_str());
            }
        }
    }

    std::string sig_path = raw_path + "/Signal";
    std::vector<i16> int_data; 
    file.read(sig_path, int_data);

    chunk_count_ = (int_data.size() / PRMS.chunk_len()) + (int_data.size() % PRMS.chunk_len() != 0);

    if (chunk_count_ > PRMS.max_chunks) {
        chunk_count_ = PRMS.max_chunks;
        int_data.resize(chunk_count_ * PRMS.chunk_len());
    }

    //full_signal_.reserve(int_data.size());

    //full_signal_.assign(int_data.begin(), int_data.end());
    //for (u16 raw : int_data) {
    for (u64 i = 0; i < samp_st; i++) {
        full_signal_.push_back(60); //DUMB
    }
    for (u64 i = samp_st; i < int_data.size(); i++) {
		float calibrated = (cal_range * int_data[i] / cal_digit) + cal_offset;
        full_signal_.push_back(calibrated);
    }

    set_raw_len(full_signal_.size());
}


ReadBuffer::ReadBuffer(Chunk &first_chunk) 
    : //source_(Source::LIVE),
      channel_idx_(first_chunk.get_channel_idx()),
      id_(first_chunk.get_id()),
      number_(first_chunk.get_number()),
      start_sample_(first_chunk.get_start()),
      chunk_count_(1),
      chunk_processed_(false) {
    set_raw_len(first_chunk.size());
    first_chunk.pop(chunk_);
}

void ReadBuffer::set_raw_len(u64 raw_len) {
    raw_len_ = raw_len;
}

bool ReadBuffer::add_chunk(Chunk &c) {
    if (!chunk_processed_ || 
        channel_idx_ != c.get_channel_idx() || 
        number_ != c.get_number()) return false;

    chunk_processed_ = false;

    chunk_count_++;
    set_raw_len(raw_len_+c.size());
    c.pop(chunk_);

    return true;
}

bool ReadBuffer::empty() const {
    return full_signal_.empty() && chunk_.empty();
}

u16 ReadBuffer::get_channel() const {
    return channel_idx_+1;
}

u16 ReadBuffer::get_channel_idx() const {
    return channel_idx_;
}

bool ReadBuffer::chunks_maxed() const {
    return chunk_count_ >= PRMS.max_chunks;
}

u32 ReadBuffer::chunk_count() const {
    return chunk_count_; //full_signal_.size() / PRMS.chunk_len();
}

Chunk ReadBuffer::get_chunk(u32 i) const {
    u32 st = i * PRMS.chunk_len(),
        ln = PRMS.chunk_len();

    if (st > full_signal_.size()) { //return Chunk();
        st = full_signal_.size();
    } 
    
    if (st+ln > full_signal_.size()) {
        ln = full_signal_.size() - st;
    }

    return Chunk(id_, get_channel(), number_, start_sample_+st, 
                 full_signal_, st, ln);
}


u32 ReadBuffer::get_chunks(std::vector<Chunk> &chunk_queue, bool real_start, u32 offs) const {
    u32 count = 0;
    u16 l = PRMS.chunk_len();

    float start = real_start ? start_sample_ : 0;

    for (u32 i = offs; i+l <= full_signal_.size() && count < PRMS.max_chunks; i += l) {
        chunk_queue.emplace_back(id_, get_channel(), number_, 
                                 start+i, full_signal_, i, l);
        count++;
    }
    return count;
}

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2) {
    return r1.start_sample_ < r2.start_sample_;
}

u64 ReadBuffer::get_duration() const {
    return raw_len_;
}

u64 ReadBuffer::get_start() const {
    return start_sample_;
}

u64 ReadBuffer::get_end() const {
    return start_sample_ + raw_len_;
}
