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

ReadParams ReadBuffer::PRMS;

const std::string Paf::PAF_TAGS[] = {
    "mt", //MAP_TIME
    "wt", //WAIT_TIME
    "qt", //QUEUE_TIME
    "rt", //RECEIVE_TIME
    "ch", //CHANNEL
    "ej", //UNBLOCK
    "st", //START_TIME
    "mx", //IN_SCAN
    "tr", //TOP_RATIO
    "mr", //MEAN_RATIO
    "en", //ENDED
    "kp", //KEEP
    "dl"  //DELAY
};

Paf::Paf() 
    : is_mapped_(false),
      ended_(false),
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

Paf::Paf(const std::string &rd_name, u16 channel, u64 start_sample)
    : is_mapped_(false),
      ended_(false),
      rd_name_(rd_name),
      rf_name_(""),
      rd_st_(0),
      rd_en_(0),
      rd_len_(0),
      rf_st_(0),
      rf_en_(0),
      rf_len_(0),
      fwd_(false),
      matches_(0) {
    
    set_int(Tag::CHANNEL, channel);
    set_int(Tag::READ_START, start_sample);
}

bool Paf::is_mapped() const {
    return is_mapped_;
}

bool Paf::is_ended() const {
    return ended_;
}

void Paf::print_paf() const {
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
        std::cout << std::fixed << "\t" << PAF_TAGS[t.first] << ":i:" << t.second;
    }
    for (auto t : float_tags_) { 
        std::cout << std::fixed << "\t" << PAF_TAGS[t.first] << ":f:" << t.second;
    }
    for (auto t : str_tags_) { 
        std::cout << "\t" << PAF_TAGS[t.first] << ":Z:" << t.second;
    }

    std::cout << "\n";
}

void Paf::set_read_len(u64 rd_len) {
    rd_len_ = rd_len;
}

void Paf::set_ended() {
    ended_ = true;
    //set_int(Tag::ENDED, 1);
}

void Paf::set_mapped(u64 rd_st, u64 rd_en,
                          std::string rf_name,
                          u64 rf_st, u64 rf_en, u64 rf_len,
                          bool fwd, u16 matches) {
    is_mapped_ = true;
    rd_st_ = rd_st;
    rd_en_ = rd_en;
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
    std::swap(loc_, r.loc_);
}

void ReadBuffer::clear() {
    raw_len_ = 0;
    full_signal_.clear();
    chunk_.clear();
    chunk_count_ = 0;
    loc_ = Paf();
}

ReadBuffer::ReadBuffer(const hdf5_tools::File &file, 
                       const std::string &raw_path, 
                       const std::string &ch_path) {

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
        } //else if (a.first == "sampling_rate") {
            //PARAMS.set_sample_rate(atof(a.second.c_str()));
        //}
    }


    set_calibration(get_channel(), offset, range, digitisation);


    std::string sig_path = raw_path + "/Signal";
    std::vector<i16> int_data; 
    file.read(sig_path, int_data);

    chunk_count_ = (int_data.size() / PRMS.chunk_len()) + (int_data.size() % PRMS.chunk_len() != 0);

    if (chunk_count_ > PRMS.max_chunks) {
        chunk_count_ = PRMS.max_chunks;
        int_data.resize(chunk_count_ * PRMS.chunk_len());
    }

    full_signal_ = calibrate(get_channel(), int_data);

    loc_ = Paf(id_, get_channel(), start_sample_);
    set_raw_len(full_signal_.size());
}


ReadBuffer::ReadBuffer(Chunk &first_chunk) 
    : //source_(Source::LIVE),
      channel_idx_(first_chunk.get_channel_idx()),
      id_(first_chunk.get_id()),
      number_(first_chunk.get_number()),
      start_sample_(first_chunk.get_start()),
      chunk_count_(1),
      chunk_processed_(false),
      loc_(id_, channel_idx_+1, start_sample_) {
    //loc_.set_int(Paf::Tag::RECEIVE_TIME, PARAMS.get_time());//TODO: FIX
    set_raw_len(first_chunk.size());
    first_chunk.pop(chunk_);
}

void ReadBuffer::set_raw_len(u64 raw_len) {
    raw_len_ = raw_len;
    loc_.set_read_len(raw_len_ * PRMS.bp_per_samp());
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

void ReadBuffer::calibrate(u16 ch, std::vector<float> samples) {
    for (float &s : samples) s = calibrate(ch, s);
}

std::vector<float> ReadBuffer::calibrate(u16 ch, std::vector<i16> samples) {
    std::vector<float> cal(samples.size());
    for (u32 i = 0; i < samples.size(); i++) {
        cal[i] = calibrate(ch, samples[i]);
    }
    return cal;
}

float ReadBuffer::calibrate(u16 ch, float sample) {
    return (sample + PRMS.calib_offsets[ch-1]) * PRMS.calib_coefs[ch-1];
    //return sample;
}

void ReadBuffer::set_calibration(const std::vector<float> &offsets, 
                                 const std::vector<float> &pa_ranges,
                                 float digitisation) {
    if (digitisation > 0) PRMS.calib_digitisation = digitisation;
    PRMS.calib_offsets = offsets;
    PRMS.calib_coefs.reserve(pa_ranges.size());
    for (float p : pa_ranges) PRMS.calib_coefs.push_back(p / digitisation);
}

void ReadBuffer::set_calibration(u16 channel,
                             float offsets, 
                             float pa_ranges,
                             float digitisation) {
    if (digitisation > 0) PRMS.calib_digitisation = digitisation;

    u16 c = channel-1;

    //TODO: will normalization take care of this?
    if (PRMS.calib_offsets.size() <= c) {
        PRMS.calib_offsets.resize(channel);
        PRMS.calib_coefs.resize(channel);
    }

    PRMS.calib_offsets[c] = offsets;
    PRMS.calib_coefs[c] = pa_ranges / digitisation;
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
