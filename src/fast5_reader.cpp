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

Chunk::Chunk() 
    : id(""),
      number(0),
      chunk_start_sample(0),
      raw_data() {}

Chunk::Chunk(const std::string &_id, u32 _number, u64 _chunk_start_sample, 
             const std::vector<float> &_raw_data, u32 raw_st=0, u32 raw_len=4000) 
    : id(_id),
      number(_number),
      chunk_start_sample(_chunk_start_sample) {
    if (raw_st + raw_len > _raw_data.size()) raw_len = _raw_data.size() - raw_st;
    raw_data.resize(raw_len);
    for (u32 i = 0; i < raw_len; i++) raw_data[i] = _raw_data[raw_st+i];
}

Chunk::Chunk(const Chunk &c) 
    : id(c.id),
      number(c.number),
      chunk_start_sample(c.chunk_start_sample),
      raw_data(c.raw_data) {}

u64 Chunk::get_end() {
    return chunk_start_sample + raw_data.size();
}

void Chunk::swap(Chunk &c) {
    std::swap(id, c.id);
    std::swap(number, c.number);
    std::swap(chunk_start_sample, c.chunk_start_sample);
    raw_data.swap(c.raw_data);
}

bool Chunk::empty() const {
    return raw_data.empty();
}

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
            channel = atoi(file.get_channel_id_params().channel_number.c_str()) - 1;
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

u32 Fast5Read::get_chunks(std::deque<Chunk> &chunk_queue, u16 max_length) {
    u32 count = 0;
    for (u32 i = 0; i < raw_data.size(); i += max_length) {
        chunk_queue.emplace_back(id, number, start_sample+i, raw_data, i, max_length);
        count++;
    }
    return count;
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
     speed_(speed),
     chshifts_(num_chs, 0),
     chunks_(num_chs) {
    timer_.reset();
    tshift_ = -1;
    is_running = true;
}

void ChunkSim::add_files(const std::vector<std::string> &fnames) {
    u32 i;
    for (i = 0; i < fnames.size() && num_loaded_ < max_loaded_; i++) {
        Fast5Read r(fnames[i]);
        r.get_chunks(chunks_[r.channel], chunk_len_);
        if (r.start_sample < tshift_) { 
            tshift_ = r.start_sample;
        }
    }
    for (; i < fnames.size(); i++) fast5_names_.push_back(fnames[i]); 
}

std::vector<ChChunk> ChunkSim::get_read_chunks() {
    u64 time = (timer_.get() * speed_) + tshift_;
    
    std::vector<ChChunk> ret;
    is_running = false;

    for (u16 c = 0; c < chunks_.size(); c++) {
        if (chunks_[c].empty()) continue;
        is_running = true; 

        //Find first chunk that ends after current time
        //Ideally will be second chunk, unless we missed some
        u16 i = 0;
        for (; i < chunks_[c].size() && 
               chunks_[c][i].chunk_start_sample+chunk_len_ < time+chshifts_[c]; i++);

        //Skip if first chunk, otherwise add previous chunk
        if (i-- == 0) {
            continue; 
        }// else if (i != 0) {
        //    std::cout << "# skipped " << i << " " << c << " "
        //              << chunks_[c][0].chunk_start_sample << "-"
        //              << chunks_[c][0].get_end() << " "
        //              << chunks_[c][1].chunk_start_sample << "-"
        //              << chunks_[c][1].get_end() << "\n";
        //}
        ret.emplace_back(c, chunks_[c][i]);

        //Remove all finished chunks
        for (; i < chunks_[c].size(); i--) {
            chunks_[c].pop_front();
        }
    }

    //if (ret.size() > 0) std::cout << time << "\t";

    return ret;
}

void ChunkSim::stop_receiving_read(u16 channel, u32 number) {
    while (!chunks_[channel].empty() && chunks_[channel][0].number == number) {
        chunks_[channel].pop_front();
    }
}

void ChunkSim::unblock(u16 channel, u32 number) {
    u64 t0 = chunks_[channel].front().chunk_start_sample;
    while (!chunks_[channel].empty() && chunks_[channel][0].number == number) {
        chunks_[channel].pop_front();
    }
    chshifts_[channel] += chunks_[channel].front().chunk_start_sample - t0;
}

