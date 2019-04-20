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

#include "fast5_reader.hpp"


u32 load_multi_fast5(const std::string &fname, std::vector<Fast5Read> &list) {
    std::string root = "/";
    hdf5_tools::File file(fname);
    std::vector<std::string> ids = file.list_group(root);
    u32 nloaded = 0;

    for (auto read : ids) {
        std::vector<int> samples;
        std::string raw_path = root + read + "/Raw",
                    sig_path = raw_path + "/Signal",
                    ch_path = root + read + "/channel_id";

        u16 channel;
        std::string id;
        u32 number;
        u64 start_sample;
        std::vector<int> int_data; 


        for (auto a : file.get_attr_map(raw_path)) {
            if (a.first == "read_id") {
                id = a.second;
            } else if (a.first == "read_number") {
                number = atoi(a.second.c_str());
            } else if (a.first == "start_time") {
                start_sample = atoi(a.second.c_str());
            }
        }

        file.read(sig_path, int_data);
        std::vector<float> raw_data(int_data.size()); 

        float digitisation = 0, range = 0, offset = 0, sampling_rate = 0;
        for (auto a : file.get_attr_map(ch_path)) {
            if (a.first == "channel_number") {
                channel = atoi(a.second.c_str()) - 1;
            } else if (a.first == "digitisation") {
                digitisation = atof(a.second.c_str());
            } else if (a.first == "range") {
                range = atof(a.second.c_str());
            } else if (a.first == "offset") {
                offset = atof(a.second.c_str());
            } else if (a.first == "sampling_rate") {
                sampling_rate = atof(a.second.c_str());
            }
        }

        if (Fast5Read::sampling_rate == 0) {
            Fast5Read::sampling_rate = sampling_rate;
        } else if (Fast5Read::sampling_rate != sampling_rate) {
            std::cerr << "Error: different sampling rates between reads, which is weird\n";
        }


        for (u32 i = 0; i < int_data.size(); i++) {
            raw_data[i] = ((float) int_data[i] + offset) * range / digitisation;
        }

        list.emplace_back(id, channel, number, start_sample, raw_data);
        nloaded++;
    }

    return nloaded;
}

float Fast5Read::sampling_rate = 0;
Fast5Read::Fast5Read(const std::string &_id,
                     u16 _channel, u32 _number, u64 _start_sample,
                     const std::vector<float> _raw_data) 
    : id(_id),
      channel(_channel),
      number(_number),
      start_sample(_start_sample),
      raw_data(_raw_data) {}

Fast5Read::Fast5Read() :
    id(""),
    channel(0),
    number(0),
    start_sample(0),
    raw_data(),
    i(0) {}

Fast5Read::Fast5Read(const std::string &filename) {

    if (filename.empty()) {
        id = "";
        channel = number = start_sample = 0;
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
            
            raw_data = file.get_raw_samples();
            
            auto p = file.get_raw_samples_params();
            id = p.read_id;
            number = p.read_number;
            start_sample = p.start_time;

            auto c = file.get_channel_id_params();
            channel = atoi(c.channel_number.c_str()) - 1;

            if (Fast5Read::sampling_rate == 0) {
                Fast5Read::sampling_rate = c.sampling_rate;
            } else if (Fast5Read::sampling_rate != c.sampling_rate) {
                std::cerr << "Error: different sampling rates between reads, which is weird\n";
            }
        }

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }
}

void Fast5Read::swap(Fast5Read &r) {
    raw_data.swap(r.raw_data);
    id.swap(r.id);
    std::swap(i, r.i);
    std::swap(channel, r.channel);
    std::swap(start_sample, r.start_sample);
    std::swap(number, r.number);
}

float Fast5Read::next_sig() {
    return raw_data[i++];
}

bool Fast5Read::empty() const {
    return i >= raw_data.size();
}

u32 Fast5Read::get_chunks(std::deque<Chunk> &chunk_queue, u16 max_length) const {
    u32 count = 0;
    for (u32 i = 0; i < raw_data.size(); i += max_length) {
        chunk_queue.emplace_back(id, channel, number, start_sample+i, raw_data, i, max_length);
        count++;
    }
    return count;
}

bool operator< (const Fast5Read &r1, const Fast5Read &r2) {
    return r1.start_sample < r2.start_sample;
}

ChunkSim::ChunkSim(u32 max_loaded, 
                   u32 num_chs, 
                   u16 chunk_len, 
                   float speed,
                   const std::vector<std::string> &fnames) 
    : ChunkSim(max_loaded, num_chs, chunk_len, speed) {
    add_files(fnames);
    timer_.reset();
}

ChunkSim::ChunkSim(u32 max_loaded, u32 num_chs, u16 chunk_len, float speed)
   : max_loaded_(max_loaded),
     num_loaded_(0),
     chunk_len_(chunk_len),
     speed_(speed / 1000),
     chshifts_(num_chs, 0),
     chunks_(num_chs) {
    tshift_ = -1;
    is_running_ = false;
}

void ChunkSim::add_reads(const std::vector<Fast5Read> &reads) {
    for (const Fast5Read &r : reads) {
        r.get_chunks(chunks_[r.channel], chunk_len_);
        if (r.start_sample < tshift_) { 
            tshift_ = r.start_sample;
        }
    }
}

void ChunkSim::add_files(const std::vector<std::string> &fnames) {
    std::vector<Fast5Read> reads;
    for (u32 i = 0; i < fnames.size(); i++) {
        load_multi_fast5(fnames[i], reads);
    }

    //TODO: sort each file? Note sure if I trust the order between files
    std::sort(reads.begin(), reads.end());
    while (reads.size() > max_loaded_) reads.pop_back();

    add_reads(reads);
}

//void ChunkSim::add_files(const std::vector<std::string> &fnames) {
//    u32 i;
//    for (i = 0; i < fnames.size() && num_loaded_ < max_loaded_; i++) {
//        Fast5Read r(fnames[i]);
//        r.get_chunks(chunks_[r.channel], chunk_len_);
//        if (r.start_sample < tshift_) { 
//            tshift_ = r.start_sample;
//        }
//    }
//    for (; i < fnames.size(); i++) fast5_names_.push_back(fnames[i]); 
//}

std::vector<Chunk> ChunkSim::get_read_chunks() {
    u64 time = (timer_.get() * speed_ * Fast5Read::sampling_rate) + tshift_;
    
    std::vector<Chunk> ret;
    is_running_ = false;

    for (u16 c = 0; c < chunks_.size(); c++) {
        if (chunks_[c].empty()) continue;
        is_running_ = true; 

        //Find first chunk that ends after current time
        //Ideally will be second chunk, unless we missed some
        u16 i = 0;
        for (; i < chunks_[c].size() && 
               chunks_[c][i].get_start()+chunk_len_ < time+chshifts_[c]; i++);

        //Skip if first chunk, otherwise add previous chunk
        if (i-- == 0) {
            continue; 
        }

        ret.push_back(chunks_[c][i]);
        //TODO: try code below
        //ret.push_back(Chunk());
        //ret.back().swap(chunks_[c][i]);

        //Remove all finished chunks
        for (; i < chunks_[c].size(); i--) {
            chunks_[c].pop_front();
        }
    }

    return ret;
}

void ChunkSim::stop_receiving_read(u16 channel, u32 number) {
    while (!chunks_[channel].empty() && 
            chunks_[channel][0].get_number() == number) {
        chunks_[channel].pop_front();
    }
}

void ChunkSim::unblock(u16 channel, u32 number) {
    u64 t0 = chunks_[channel].front().get_start();
    while (!chunks_[channel].empty() && 
            chunks_[channel][0].get_number() == number) {
        chunks_[channel].pop_front();
    }
    chshifts_[channel] += chunks_[channel].front().get_start() - t0;
}

void ChunkSim::set_time(ReadLoc &loc) {
    u16 c = loc.get_channel();
    u64 time = (timer_.get() * speed_ * Fast5Read::sampling_rate) + tshift_ + chshifts_[c];
    loc.set_time(time);
}

void ChunkSim::start() {
    is_running_ = true;
    timer_.reset();
}

bool ChunkSim::is_running() {
    return is_running_;
}
