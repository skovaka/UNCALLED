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

#include "paf.hpp"

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
    "dl", //DELAY
    "sc", //SEED_CLUSTER
    "ce"  //CONFIDENT_EVENT
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


