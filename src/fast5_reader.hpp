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

#ifndef FAST5_READER_HPP
#define FAST5_READER_HPP

#include <deque>
#include "fast5.hpp"
#include "util.hpp"
#include "timer.hpp"

class Chunk {
    public:
    Chunk(const std::string &_id, u32 _number, u64 _chunk_start_sample, 
          const std::vector<float> &_raw_data, u32 raw_st, u32 raw_len);
    Chunk(const Chunk &c);

    u64 get_end();

    std::string id;
    u32 number;
    u64 chunk_start_sample;
    std::vector<float> raw_data;
    //std::vector<u32> chunk_classifications;
    //float median_before, median;

};
using ChChunk = std::pair<u16, Chunk>;

class Fast5Read {
    public:
    Fast5Read();
    Fast5Read(const std::string &filename);

    //std::vector<Chunk> get_chunks(u16 max_length);
    u32 get_chunks(std::deque<Chunk> &chunk_queue, u16 max_length);
    void swap(Fast5Read &r);
    float next_sig();
    bool empty() const;

    std::string id;
    u16 channel;
    u32 number;
    u64 start_sample;
    std::vector<float> raw_data;
    u32 i;
    //float median_before;
};


class ChunkSim {
    public:
    ChunkSim(u32 max_loaded, u32 num_chs, u16 chunk_len, const std::vector<std::string> &fnames);
    ChunkSim(u32 max_loaded, u32 num_chs, u16 chunk_len);
    
    void add_files(const std::vector<std::string> &fnames);

    std::vector<ChChunk> get_read_chunks();
    void stop_receiving_read(u16 channel, u32 number);

    private:
    u32 max_loaded_, num_loaded_;
    u16 chunk_len_;
    Timer timer_;
    std::deque<std::string> fast5_names_;
    std::vector< std::deque<Chunk> > chunks_;
};

#endif
