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

#ifndef READ_BUFFER_HPP
#define READ_BUFFER_HPP

#include <string>
#include <vector>
#include "fast5/hdf5_tools.hpp"
#include "util.hpp"
#include "chunk.hpp"

class Paf {
    public:
    enum Tag {MAP_TIME, CHUNKS, UNBLOCK, KEEP};

    Paf();
    Paf(std::string rd_name_, u64 rd_len = 0);

    bool is_mapped() const;
    void print_paf() const;
    void set_read_len(u64 rd_len);
    void set_mapped(u64 rd_st, u64 rd_en, 
                    std::string rf_name,
                    u64 rf_st, u64 rf_en, u64 rf_len,
                    bool fwd, u16 matches);
    void set_unmapped();

    void set_int(Tag t, int v);
    void set_float(Tag t, float v);
    void set_str(Tag t, std::string v);

    private:
    static const std::string PAF_TAGS[];

    bool is_mapped_;
    std::string rd_name_, rf_name_;
    u64 rd_st_, rd_en_, rd_len_,
        rf_st_, rf_en_, rf_len_;
    bool fwd_;
    u16 matches_;

    std::vector< std::pair<Tag, int> > int_tags_;
    std::vector< std::pair<Tag, float> > float_tags_;
    std::vector< std::pair<Tag, std::string> > str_tags_;
};

class ReadBuffer {
    public:

    enum Source {MULTI, SINGLE, BULK, LIVE};


    ReadBuffer();
    ReadBuffer(const ReadBuffer &read);
    ReadBuffer(const std::string &filename);
    ReadBuffer(const hdf5_tools::File &file, Source source, const std::string root="/");
    ReadBuffer(Source source, u16 channel, const std::string &id = "", 
               u32 number = 0, u64 start_sample = 0, 
               const std::vector<float> raw_data = std::vector<float>(),
               u32 raw_st = 0, u32 raw_len = 0);
    ReadBuffer(Chunk &first_chunk);

    void fast5_init(const hdf5_tools::File &file, 
                    std::string raw_path, 
                    std::string ch_path);

    bool add_chunk(Chunk &c);

    void print_paf();


    Chunk &&pop_chunk();
    void swap(ReadBuffer &r);
    void clear();
    bool empty() const;

    void set_raw_len(u64 raw_len_);

    u32 get_chunks(std::deque<Chunk> &chunk_queue, u16 max_length) const;

    u16 get_channel() const;
    u16 get_channel_idx() const;

    static void set_calibration(float digitisation, 
                                const std::vector<float> &offsets, 
                                const std::vector<float> &pa_ranges);
    static float sampling_rate;
    static std::vector<float> cal_offsets_, cal_coefs_;

    Source source_;
    u16 channel_idx_;
    std::string id_;
    u32 number_;
    u64 start_sample_, raw_len_;
    std::vector<float> full_signal_, chunk_;
    u16 num_chunks_;
    bool chunk_processed_;

    Paf loc_;

    friend bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);
};

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);

u32 load_fast5s(const std::string &fname, std::deque<ReadBuffer> &list, u32 max_load = 0);
u32 load_fast5s(std::ifstream &file, std::deque<ReadBuffer> &list, u32 max_load = 0);
u32 load_multi_fast5(const hdf5_tools::File &file, std::deque<ReadBuffer> &list, u32 max_load = 0);

#endif
