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
#include "params.hpp"

const std::string Paf::PAF_TAGS[] = {
    "mt", //MAP_TIME
    "nc", //CHUNKS
    "ub", //UNBLOCK
    "kp"  //KEEP
};

Paf::Paf() 
    : is_mapped_(false),
      rd_name_(""),
      rf_name_(""),
      rd_st_(0),
      rd_en_(0),
      rd_len_(0),
      rf_st_(0),
      rf_en_(0),
      rf_len_(0),
      fwd_(false),
      matches_(0) {}

Paf::Paf(std::string rd_name, u64 rd_len) 
    : is_mapped_(false),
      rd_name_(rd_name),
      rf_name_(""),
      rd_st_(0),
      rd_en_(0),
      rd_len_(rd_len),
      rf_st_(0),
      rf_en_(0),
      rf_len_(0),
      fwd_(false),
      matches_(0) {}

bool Paf::is_mapped() const {
    return is_mapped_;
}

void Paf::print() const {
    std::cout << rd_name_ << "\t"
       << rd_len_ << "\t";
    if (is_mapped_) {
       std::cout 
           << rd_st_ << "\t"
           << rd_en_ << "\t"
           << (fwd_ ? '+' : '-') << "\t"
           << rf_name_ << "\t"
           << rf_len_ << "\t"
           << rf_st_ << "\t"
           << rf_en_ << "\t"
           << matches_ << "\t"
           << (rf_en_ - rf_st_ + 1) << "\t"
           << 255;
    } else {
        std::cout << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "255";
    }

    for (auto t : int_tags_) { 
        std::cout << "\t" << PAF_TAGS[t.first] << ":i:" << t.second;
    }
    for (auto t : float_tags_) { 
        std::cout << "\t" << PAF_TAGS[t.first] << ":f:" << t.second;
    }
    for (auto t : str_tags_) { 
        std::cout << "\t" << PAF_TAGS[t.first] << ":Z:" << t.second;
    }

    std::cout << std::endl;
}

void Paf::set_read_len(u64 rd_len) {
    rd_len_ = rd_len;
}

void Paf::set_mapped(u64 rd_st, u64 rd_en, u64 rd_len,
                          std::string rf_name,
                          u64 rf_st, u64 rf_en, u64 rf_len,
                          bool fwd, u16 matches) {
    is_mapped_ = true;
    rd_st_ = rd_st;
    rd_en_ = rd_en;
    rd_len_ = rd_len;
    rf_name_ = rf_name;
    rf_st_ = rf_st;
    rf_en_ = rf_en;
    rf_len_ = rf_len;
    fwd_ = fwd;
    matches_ = matches;
}

void Paf::set_int(Tag t, int v) {
    int_tags_.emplace_back(t, v);
}

void Paf::set_float(Tag t, float v) {
    float_tags_.emplace_back(t, v);
}

void Paf::set_str(Tag t, std::string v) {
    str_tags_.emplace_back(t, v);
}

bool is_multi_fast5(const hdf5_tools::File &file) {
    for (const std::string &s : file.list_group("/")) {
        if (s == "Raw") return false;
    }
    return true;
}

u32 load_multi_fast5(const hdf5_tools::File &file, std::vector<ReadBuffer> &list, u32 max_load) { 
    u32 i = 0;
    for (const std::string &read : file.list_group("/")) {
        list.emplace_back(file, ReadBuffer::Source::MULTI, "/" + read);
        i++;
        if (max_load > 0 && list.size() >= max_load) break;
    }
    return i;
}

u32 load_fast5s(const std::string &fname, 
                std::vector<ReadBuffer> &list, 
                u32 max_load) {
    std::string root = "/";

    std::ifstream reads_file(fname);
    if (!reads_file) {
        std::cerr << "Error: couldn't open '" << fname << "'\n";
        return 1;
    }

    u32 i = 0;
    std::string read_fname;
    while (getline(reads_file, read_fname) && 
            (max_load == 0 || list.size() < max_load)) {
        if (read_fname[0] == '#') continue;
        hdf5_tools::File file(read_fname);
        
        if (is_multi_fast5(file)) {
            i += load_multi_fast5(file, list, max_load);
        } else {
            list.emplace_back(file, ReadBuffer::Source::SINGLE);
            i++;
        }
    }

    return i;
}

ReadBuffer::ReadBuffer() {
    
}

ReadBuffer::ReadBuffer(const ReadBuffer &r) 
    : source_          (r.source_),
      channel_idx_     (r.channel_idx_),
      id_              (r.id_),
      number_          (r.number_),
      start_sample_    (r.start_sample_),
      raw_len_         (r.raw_len_),
      full_signal_     (r.full_signal_),
      chunk_           (r.chunk_),
      num_chunks_      (r.num_chunks_),
      chunk_processed_ (r.chunk_processed_),
      loc_             (r.loc_) {}

ReadBuffer::ReadBuffer(const hdf5_tools::File &file, Source source, 
                       const std::string root) {
    source_ = source;
    std::string raw_path, ch_path;
                 
    switch (source) {
    case Source::MULTI:
        raw_path = root + "/Raw";
        ch_path = root + "/channel_id";
        fast5_init(file, raw_path, ch_path);
    break;
    case Source::SINGLE:
        raw_path = root + "Raw/Reads";
        ch_path = root + "UniqueGlobalKey/channel_id";
        raw_path += "/" + file.list_group(raw_path).front();
        fast5_init(file, raw_path, ch_path);
    break;

    case Source::BULK:
        std::cerr << "Bulk not supported yet\n";
    break;
    
    default:
        std::cerr << "Error: invalid fast5 source\n";
    break;
    }

    loc_ = Paf(id_);
}

ReadBuffer::ReadBuffer(Source source, u16 channel, const std::string &id, 
                       u32 number, u64 start_sample, 
                       const std::vector<float> raw_data,
                       u32 raw_st, u32 raw_len) 
        : source_(source),
          channel_idx_(channel-1),
          id_(id),
          number_(number),
          start_sample_(start_sample),
          loc_(id) {
    if (raw_data.empty()) raw_len_ = 0;
    else {
        if (raw_len == 0) raw_len = raw_data.size() - raw_st;
        full_signal_ = std::vector<float>(&raw_data[raw_st], 
                                          &raw_data[raw_st+raw_len]);
        raw_len_ = raw_data.size();
    }
}

ReadBuffer::ReadBuffer(Chunk &first_chunk) 
    : source_(Source::LIVE),
      channel_idx_(first_chunk.get_channel_idx()),
      id_(first_chunk.get_id()),
      number_(first_chunk.get_number()),
      start_sample_(first_chunk.get_start()),
      raw_len_(first_chunk.size()),
      num_chunks_(1),
      chunk_processed_(false),
      loc_(id_) {
    first_chunk.pop(chunk_);
}

void ReadBuffer::fast5_init(const hdf5_tools::File &file, 
                            std::string raw_path, 
                            std::string ch_path) {


    for (auto a : file.get_attr_map(raw_path)) {
        if (a.first == "read_id") {
            id_ = a.second;
        } else if (a.first == "read_number") {
            number_ = atoi(a.second.c_str());
        } else if (a.first == "start_time") {
            start_sample_ = atoi(a.second.c_str());
        }
    }

    float digitisation = 0, range = 0, offset = 0;//, sampling_rate = 0;
    for (auto a : file.get_attr_map(ch_path)) {
        if (a.first == "channel_number") {
            channel_idx_ = atoi(a.second.c_str()) - 1;
        } else if (a.first == "digitisation") {
            digitisation = atof(a.second.c_str());
        } else if (a.first == "range") {
            range = atof(a.second.c_str());
        } else if (a.first == "offset") {
            offset = atof(a.second.c_str());
        } else if (a.first == "sampling_rate") {
            PARAMS.set_sample_rate(atof(a.second.c_str()));
        }
    }

    PARAMS.set_calibration(get_channel(), offset, range, digitisation);

    std::string sig_path = raw_path + "/Signal";
    std::vector<i16> int_data; 
    file.read(sig_path, int_data);
    full_signal_ = PARAMS.calibrate(get_channel(), int_data);

    /*
    full_signal_.resize(int_data.size());
    for (u32 i = 0; i < int_data.size(); i++) {
        //std::cout << int_data[i] << " ";
        full_signal_[i] = ((float) int_data[i] + offset) * range / digitisation;
        //std::cout << full_signal_[i] << "\n";
    }*/
}

bool ReadBuffer::add_chunk(Chunk &c) {
    if (!chunk_processed_ || 
        channel_idx_ != c.get_channel_idx() || 
        number_ != c.get_number()) return false;
    num_chunks_++;
    raw_len_ += c.size();
    chunk_processed_ = false;
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

u32 ReadBuffer::get_chunks(std::deque<Chunk> &chunk_queue, u16 max_length) const {
    u32 count = 0; //
    for (u32 i = 0; i+max_length <= full_signal_.size(); i += max_length) {
        chunk_queue.emplace_back(id_, get_channel(), number_, 
                                 start_sample_+i, full_signal_, 
                                 i, max_length);
        count++;
    }
    return count;
}

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2) {
    return r1.start_sample_ < r2.start_sample_;
}

