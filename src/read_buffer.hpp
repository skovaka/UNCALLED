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
#include <unordered_set>
#include <fast5/hdf5_tools.hpp>
#include "util.hpp"
#include "chunk.hpp"


typedef struct {
    u16 num_channels;
    float bp_per_sec;
    float sample_rate;
    float calib_digitisation;
    float chunk_time;
    u32 max_chunks;
    std::vector<float> calib_offsets, calib_coefs;

    float bp_per_samp() {
        return bp_per_sec / sample_rate;
    }

    u16 chunk_len() {
        return (u16) (chunk_time * sample_rate);
    }
} ReadParams;

class Paf {
    public:

    enum Tag {MAP_TIME, 
              WAIT_TIME, 
              RECEIVE_TIME,
              CHANNEL, 
              EJECT, 
              READ_START, 
              IN_SCAN, 
              TOP_RATIO, 
              MEAN_RATIO,
              ENDED,
              KEEP,
              DELAY};

    Paf();
    Paf(const std::string &rd_name, u16 channel = 0, u64 start_sample = 0);

    bool is_mapped() const;
    bool is_ended() const;
    void print_paf() const;
    void set_read_len(u64 rd_len);
    void set_mapped(u64 rd_st, u64 rd_en, 
                    std::string rf_name,
                    u64 rf_st, u64 rf_en, u64 rf_len,
                    bool fwd, u16 matches);
    void set_ended();
    void set_unmapped();

    void set_int(Tag t, int v);
    void set_float(Tag t, float v);
    void set_str(Tag t, std::string v);


    private:
    static const std::string PAF_TAGS[];

    bool is_mapped_, ended_;
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
    static ReadParams PRMS;

    //enum Source {MULTI, SINGLE, BULK, LIVE};

    ReadBuffer();
    //ReadBuffer(const ReadBuffer &read);
    ReadBuffer(const std::string &filename);
    ReadBuffer(const hdf5_tools::File &file, const std::string &raw_path, const std::string &ch_path);

    //ReadBuffer(Source source, u16 channel, const std::string &id = "", 
    //           u32 number = 0, u64 start_sample = 0, 
    //           const std::vector<float> raw_data = std::vector<float>(),
    //           u32 raw_st = 0, u32 raw_len = 0);
    
    ReadBuffer(Chunk &first_chunk);

    bool empty() const;
    std::string get_id() const {return id_;}
    u64 get_start() const;
    u64 get_end() const;
    u64 get_duration() const;
    u32 size() const {return full_signal_.size();}
    u16 get_channel() const;
    const std::vector<float> &get_raw() const {return full_signal_;}

    bool add_chunk(Chunk &c);
    Chunk &&pop_chunk();
    void swap(ReadBuffer &r);
    void clear();
    void set_raw_len(u64 raw_len_);

    u32 chunk_count() const;
    bool chunks_maxed() const ;
    Chunk get_chunk(u32 i) const;

    u32 get_chunks(std::vector<Chunk> &chunk_queue, bool real_start=true, u32 offs=0) const;
    void set_channel(u16 ch) {channel_idx_ = ch-1;}
    u16 get_channel_idx() const;

    u32 get_number() const {
        return number_;
    }

    static std::vector<float> calibrate(u16 ch, std::vector<i16> samples);
    static void calibrate(u16 ch, std::vector<float> samples);
    static float calibrate(u16 ch, float sample);
    static void set_calibration(const std::vector<float> &offsets, const std::vector<float> &pa_ranges, float digitisation);
    static void set_calibration(u16 channel, float offsets, float pa_ranges, float digitisation);

    //Source source_;
    u16 channel_idx_;
    std::string id_;
    u32 number_;
    u64 start_sample_, raw_len_;
    std::vector<float> full_signal_, chunk_;
    u16 chunk_count_;
    bool chunk_processed_;

    Paf loc_;

    friend bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);
};

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);

#endif
