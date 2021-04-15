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
    fast5_list : {},
    read_filter  : {},
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
    fast5_idx_ = 0;
    total_buffered_ = 0;

    for (auto &fname : p.fast5_list) {
        add_fast5(fname);
    }
    for (auto &id : p.read_filter) {
        add_read(id);
    }
}

Fast5Reader::Fast5Reader(const std::string &fast5_list, 
                         const std::string &read_filter)
    : Fast5Reader() {

    total_buffered_ = 0;
    if (!PRMS.fast5_list.empty()) load_fast5_list(fast5_list);
    if (!PRMS.read_filter.empty()) load_read_filter(read_filter);
}

Fast5Reader::Fast5Reader(
        const std::vector<std::string> &fast5s, 
        const std::vector<std::string> &reads, 
        const Params &p) : Fast5Reader(p) {
    for (auto &fname : fast5s) {
        add_fast5(fname);
    }
    for (auto &id : reads) {
        add_read(id);
    }
}

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

bool Fast5Reader::load_read_filter(const std::string &fname) {
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
           (fast5_idx_ >= fast5_list_.size() || all_buffered());
}

bool Fast5Reader::open_fast5(u32 i) {

    if (fast5_file_.is_open()) fast5_file_.close();
    read_paths_.clear();

    //if (fast5_list_.empty()) return false;
    if (i >= fast5_list_.size()) return false;

    //fast5_file_.open(fast5_list_.front());
    //fast5_list_.pop_front();

    fast5_file_.open(fast5_list_[i]);

    fast5_fmt_ = Format::UNKNOWN;
    if (fast5_file_.group_exists(SINGLE_RAW_PATH)) {
        fast5_fmt_ = Format::SINGLE;
    } else {
        fast5_fmt_ = Format::MULTI; //TODO: add support for old multi format
    }

    return load_read_paths();
}

Fast5Read::Paths Fast5Reader::get_subpaths(const std::string &base) {
    Fast5Read::Paths paths;

    switch (fast5_fmt_) {
        case Format::SINGLE:
            paths.raw      = SINGLE_RAW_PATH + base;
            paths.channel  = SINGLE_CH_PATH;
            paths.analysis = ANALYSIS_PATH;
            break;
        case Format::MULTI:
            paths.raw      = base + MULTI_RAW_PATH;
            paths.channel  = base + MULTI_CH_PATH;
            paths.analysis = base + ANALYSIS_PATH;
            break;
        default:
            std::cerr << "Error: unrecognized fast5 format\n";
    }

    if (!PRMS.load_bc) {
        paths.analysis = "";
    }
    
    return paths;
}

bool Fast5Reader::load_read_paths() {
    std::string path;
    switch (fast5_fmt_) {
    case Format::SINGLE:

        path = SINGLE_RAW_PATH;
        for (const std::string &read : fast5_file_.list_group(path)) {

            //TODO use it like a hash table
            std::string read_id = "";
            for (auto a : fast5_file_.get_attr_map(path+"/"+read)) {
                if (a.first == "read_id") {
                    read_id = a.second;
                    break;
                }
            }
            if (read_id.empty()) {
                std::cerr << "Error: failed to find read_id\n";
                return false;
            }
            
            //TODO combine with below and put at end
            if (read_filter_.empty() || read_filter_.count(read_id) > 0) {
                read_paths_.push_back("/"+read);
            }
        }
        return true;

    //TODO put in function
    case Format::MULTI:
        path = MULTI_RAW_PATH;
        for (const std::string &read : fast5_file_.list_group("/")) {
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

        //Open fast5s until one is found which contains reads to load
        while (read_paths_.empty()) {
            if(!open_fast5(fast5_idx_++)) { 
                break;
            }
        }
        if (read_paths_.empty()) break;

        auto subpaths = get_subpaths(read_paths_.front());
        buffered_reads_.emplace_back(fast5_file_, subpaths);
        read_paths_.pop_front();


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

