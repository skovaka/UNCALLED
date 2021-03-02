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

#include <thread>
#include <chrono>
#include "fast5_reader.hpp"


const typename Fast5Reader::Params Fast5Reader::PRMS_DEF = {
    fast5_list : "",
    read_list  : "",
    max_reads  : 0,
    max_buffer : 100,
    load_bc    : false
};

const std::string 
    MULTI_RAW_PATH   = "/Raw", 
    MULTI_CH_PATH  = "/channel_id",

    SINGLE_RAW_PATH  = "/Raw/Reads",
    SINGLE_CH_PATH = "/UniqueGlobalKey/channel_id",

    ANALYSIS_PATH    = "/Analyses";

Fast5Reader::Fast5Reader() : 
    Fast5Reader(PRMS_DEF) {}

Fast5Reader::Fast5Reader(const Params &p) : PRMS(p) {
    total_buffered_ = 0;

    if (!PRMS.read_list.empty()) load_read_list(PRMS.read_list);
    if (!PRMS.fast5_list.empty()) load_fast5_list(PRMS.fast5_list);
}

Fast5Reader::Fast5Reader(const std::string &fast5_list, 
                         const std::string &read_list,
                         u32 max_reads, u32 max_buffer) 
    : PRMS({fast5_list, 
            read_list, 
            max_reads, 
            max_buffer}) {

    total_buffered_ = 0;
    if (!PRMS.fast5_list.empty()) load_fast5_list(PRMS.fast5_list);
    if (!PRMS.read_list.empty()) load_read_list(PRMS.read_list);
}

Fast5Reader::Fast5Reader(const std::vector<std::string> &fast5s, const std::vector<std::string> &reads, const Params &p) : Fast5Reader(p) {
    for (auto &fname : fast5s) {
        add_fast5(fname);
    }
    for (auto &id : reads) {
        add_read(id);
    }
}

Fast5Reader::Fast5Reader(u32 max_reads, u32 max_buffer) 
    : Fast5Reader("","",max_reads,max_buffer) {}

void Fast5Reader::add_fast5(const std::string &fast5_path) {
    fast5_list_.push_back(fast5_path);
}

bool Fast5Reader::load_fast5_list(const std::string &fname) {
    std::ifstream list_file(fname);

    if (!list_file.is_open()) {
        std::cerr << "Error: failed to open fast5 list \""
                  << fname << "\".\n";
        return false;
    }

    std::string fast5_name;
    while (getline(list_file, fast5_name)) {
        add_fast5(fast5_name);
    }

    return true;
}

bool Fast5Reader::add_read(const std::string &read_id) {
    if (PRMS.max_reads != 0 && read_filter_.size() >= PRMS.max_reads)
        return false;

    read_filter_.insert(read_id);

    return true;
}

bool Fast5Reader::load_read_list(const std::string &fname) {
    std::ifstream read_file(fname);

    if (!read_file.is_open()) {
        std::cerr << "Error: failed to open read list \""
                  << fname << "\".\n";
        return false;
    }

    std::string read_name;

    while (getline(read_file, read_name)) {
        if (!add_read(read_name)) break;
    }
    
    return true;
}

bool Fast5Reader::empty() {
    return buffered_reads_.empty() && 
           read_paths_.empty() && 
           (fast5_list_.empty() || all_buffered());
}

bool Fast5Reader::open_next() {

    read_paths_.clear();
    if (open_fast5_.is_open()) open_fast5_.close();
    if (fast5_list_.empty()) return false;

    open_fast5_.open(fast5_list_.front());
    fast5_list_.pop_front();

    open_fmt_ = Format::UNKNOWN;
    for (const std::string &s : open_fast5_.list_group("/")) {
        if (s == "Raw") {
            open_fmt_ = Format::SINGLE;
            break;
        }
    }
    if (open_fmt_ == Format::UNKNOWN) open_fmt_ = Format::MULTI; //TODO: add support for old multi format


    std::string path;
    switch (open_fmt_) {
    case Format::SINGLE:

        //TODO put in function, store attr map?
        path = SINGLE_RAW_PATH;
        for (const std::string &read : open_fast5_.list_group(path)) {
            std::string read_id = "";
            for (auto a : open_fast5_.get_attr_map(path+"/"+read)) {
                if (a.first == "read_id") {
                    read_id = a.second;
                    break;
                }
            }

            if (read_id.empty()) {
                std::cerr << "Error: failed to find read_id\n";
                return false;
            }
            
            if (read_filter_.empty() || read_filter_.count(read_id) > 0) {
                read_paths_.push_back("/"+read);
            }
        }
        return true;

    //TODO put in function
    case Format::MULTI:
        path = MULTI_RAW_PATH;
        for (const std::string &read : open_fast5_.list_group("/")) {
            std::string id = read.substr(read.find('_')+1);
            if (read_filter_.empty() || read_filter_.count(id) > 0) {
                read_paths_.push_back("/"+read);
            }
        }
        return true;
    default:
        return false;
    }

    return false; 
}

u32 Fast5Reader::fill_buffer() {
    u32 count = 0;

    //TODO: max total default to max int
    while (buffered_reads_.size() < PRMS.max_buffer) {

        if (all_buffered())  {
            read_paths_.clear();
            fast5_list_.clear();
            break;
        }

        while (read_paths_.empty()) {
            if(!open_next()) break;
        }
        if (read_paths_.empty()) break;

        std::string raw_path, ch_path, anl_path;

        switch (open_fmt_) {
            case Format::SINGLE:
                raw_path = SINGLE_RAW_PATH + read_paths_.front();
                ch_path = SINGLE_CH_PATH;
                anl_path = ANALYSIS_PATH;
                break;
            case Format::MULTI:
                raw_path = read_paths_.front() + MULTI_RAW_PATH;
                ch_path =  read_paths_.front() + MULTI_CH_PATH;
                anl_path = read_paths_.front() + ANALYSIS_PATH;
                break;
            default:
                std::cerr << "Error: unrecognized fast5 format\n";
                read_paths_.pop_front();
                return count;
        }

        if (!PRMS.load_bc) {
            anl_path = "";
        }

        read_paths_.pop_front();

        buffered_reads_.emplace_back(open_fast5_, raw_path, ch_path, anl_path);

        count++;
        total_buffered_++;
    }

    return count;
}

bool Fast5Reader::all_buffered() {
    return (PRMS.max_reads > 0 && total_buffered_ >= PRMS.max_reads) ||
           (!read_filter_.empty() && total_buffered_ >= read_filter_.size());
}

Fast5Read Fast5Reader::pop_read() {
    if (buffer_size() == 0) { 
        fill_buffer();
    }
    //TODO: swap to speed up?
    auto r = buffered_reads_.front();
    buffered_reads_.pop_front();
    return r;
}

u32 Fast5Reader::buffer_size() {
    return buffered_reads_.size();
}

