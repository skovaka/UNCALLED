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
    fast5_files : {},
    read_filter  : {},
    fast5_index : "",
    max_reads  : 0,
    max_buffer : 100,
    recursive  : false,
    load_bc    : false,
    bc_group : "000",
};

const std::string 
    MULTI_RAW_PATH   = "/Raw", 
    MULTI_CH_PATH  = "/channel_id",
    MULTI_READ_PREFIX = "read_",

    SINGLE_RAW_PATH  = "/Raw/Reads",
    SINGLE_CH_PATH = "/UniqueGlobalKey/channel_id",

    //ANALYSIS_PATH    = "",

    GUPPY_SEG_PREFIX = "/Analyses/Segmentation_",
    GUPPY_BC_PREFIX = "/Analyses/Basecall_1D_";

Fast5Reader::Fast5Reader() : 
    Fast5Reader(PRMS_DEF) {}

Fast5Reader::Fast5Reader(const Params &p) : PRMS(p) {}

u32 Fast5Reader::add_fast5(const std::string &fast5_path) {
    u32 i = fast5_files_.size();
    fast5_files_.push_back(fast5_path);
    return i;
}

bool Fast5Reader::load_fast5_files(const std::string &fname) {
    fast5_files_ = read_txt_file(fname);
    return !fast5_files_.empty();
}

bool Fast5Reader::open_fast5(u32 i) {

    if (fast5_file_.is_open()) fast5_file_.close();
    read_paths_.clear();

    //if (fast5_files_.empty()) return false;
    if (i >= fast5_files_.size()) return false;

    fast5_file_.open(fast5_files_[i]);

    fmt_ = Format::UNKNOWN;
    if (fast5_file_.group_exists(SINGLE_RAW_PATH)) {
        fmt_ = Format::SINGLE;
    } else {
        fmt_ = Format::MULTI; //TODO: add support for old multi format
    }

    return true;
}

Fast5Read::Paths Fast5Reader::get_subpaths(const std::string &path) {
    Fast5Read::Paths subpaths;

    std::string analysis_prefix;

    switch (fmt_) {
        case Format::SINGLE:
            subpaths.raw     = path;
            subpaths.channel = SINGLE_CH_PATH;
            analysis_prefix  = "";
            break;
        case Format::MULTI:
            subpaths.raw     = path + MULTI_RAW_PATH;
            subpaths.channel = path + MULTI_CH_PATH;
            analysis_prefix  = path;
            break;
        default:
            std::cerr << "Error: unrecognized fast5 format\n";
    }

    if (PRMS.load_bc) {
        subpaths.basecall = analysis_prefix + GUPPY_BC_PREFIX + PRMS.bc_group;
        subpaths.segmentation = analysis_prefix +GUPPY_SEG_PREFIX + PRMS.bc_group;
    } else {
        subpaths.basecall = subpaths.segmentation = "";
    }
    
    return subpaths;
}

Fast5Iter::Fast5Iter() : Fast5Reader() {
    total_buffered_ = 0;
    fast5_idx_ = 0;
}

Fast5Iter::Fast5Iter(const Params &p) : Fast5Reader(p) {
    total_buffered_ = 0;
    fast5_idx_ = 0;

    for (auto &fname : p.fast5_files) {
        add_fast5(fname);
    }
    for (auto &id : p.read_filter) {
        add_read(id);
    }
}

Fast5Iter::Fast5Iter(const std::string &fast5_files, 
                         const std::string &read_filter)
    : Fast5Iter() {

    total_buffered_ = 0;
    if (!PRMS.fast5_files.empty()) load_fast5_files(fast5_files);
    if (!PRMS.read_filter.empty()) load_read_filter(read_filter);
}

Fast5Iter::Fast5Iter(
        const std::vector<std::string> &fast5s, 
        const std::vector<std::string> &reads, 
        const Params &p) : Fast5Iter(p) {
    for (auto &fname : fast5s) {
        add_fast5(fname);
    }
    for (auto &id : reads) {
        add_read(id);
    }
}


bool Fast5Iter::add_read(const std::string &read_id) {
    if (PRMS.max_reads != 0 && read_filter_.size() >= PRMS.max_reads)
        return false;

    read_filter_.insert(read_id);

    return true;
}

bool Fast5Iter::load_read_filter(const std::string &fname) {
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

bool Fast5Iter::empty() {
    if (buffer_size() > 0) { 
        return false;
    }
    return fill_buffer() == 0;

    //return buffered_reads_.empty() && 
    //       read_paths_.empty() && 
    //       (fast5_idx_ >= fast5_files_.size() || all_buffered());
}

std::string Fast5Reader::get_single_raw_path() {
    return SINGLE_RAW_PATH + "/" + fast5_file_.list_group(SINGLE_RAW_PATH)[0];
}

std::string Fast5Reader::get_read_id(const std::string &raw_path) {
    return fast5_file_.get_attr_map(raw_path)["read_id"];
}

//bool Fast5Iter::load_read_paths() {
bool Fast5Iter::open_fast5(u32 i) {
    if (!Fast5Reader::open_fast5(i)) {
        return false;
    }

    std::string read_id, path;
    switch (fmt_) {
    case Format::SINGLE:
        path = get_single_raw_path();
        read_id = get_read_id(path);

        if (read_filter_.empty() || read_filter_.count(read_id) > 0) {
            read_paths_.push_back(path);
        }
        return true;

    //TODO put in function
    case Format::MULTI:
        for (const std::string &read_path : fast5_file_.list_group("/")) {
            
            read_id = read_path.substr(MULTI_READ_PREFIX.size());
            if (read_filter_.empty() || read_filter_.count(read_id) > 0) {
                read_paths_.push_back("/"+read_path);
            }
        }
        return true;

    default:
        return false;
    }

    return false; 
}

u32 Fast5Iter::fill_buffer() {
    u32 count = 0;

    //TODO: max total default to max int
    while (buffered_reads_.size() < PRMS.max_buffer) {

        if (all_buffered())  {
            read_paths_.clear();
            fast5_files_.clear();
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

        try {
            buffered_reads_.emplace_back(fast5_file_, subpaths);
        } catch (const std::exception e) {
            std::cerr << "Failed parse read in file \""
                      << read_paths_.front() << "\" at \""
                      << subpaths.raw << "\"\n";
        }

        read_paths_.pop_front();

        count++;
        total_buffered_++;
    }

    return count;
}

bool Fast5Iter::all_buffered() {
    return (PRMS.max_reads > 0 && total_buffered_ >= PRMS.max_reads) ||
           (!read_filter_.empty() && total_buffered_ >= read_filter_.size());
}

Fast5Read Fast5Iter::next_read() {
    if (empty()) {
        throw std::runtime_error("Fast5Iter is empty");
    }

    auto r = buffered_reads_.front();
    buffered_reads_.pop_front();
    return r;
}

u32 Fast5Iter::buffer_size() {
    return buffered_reads_.size();
}

Fast5Dict::Fast5Dict() : Fast5Dict(PRMS_DEF) {}

Fast5Dict::Fast5Dict(const Params &p) : Fast5Reader(p) {
    fast5_idx_ = -1;
    if (!PRMS.fast5_index.empty()) {
        load_index(PRMS.fast5_index);
    }
}

Fast5Dict::Fast5Dict(std::vector<std::string> file_paths, std::vector<std::string> read_ids, std::vector<std::string> file_names, const Params &p) : Fast5Dict(p) {
    for (auto &path : file_paths) {
        auto i = add_fast5(path);
        auto basename = path.substr(path.rfind('/')+1);
        filename_paths_[basename] = i;
    }   

    if (read_ids.size() != file_names.size()) {
        throw std::runtime_error("read_ids must be same length as file_names");
    }

    for (size_t i = 0; i < read_ids.size(); i++) {
        add_read(read_ids[i], filename_paths_[file_names[i]]);
    }

}

Fast5Dict::Fast5Dict(Fast5ReadMap fast5_map, const Params &p) : Fast5Dict(p) {
    for (auto fast5_reads : fast5_map) {
        auto i = add_fast5(fast5_reads.first);
        for (auto &read : fast5_reads.second) {
            add_read(read, i);
        }
    }
}

bool Fast5Dict::load_index(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed to open \"" << filename << "\".\n";
        return false;
    }

    std::unordered_map<std::string, size_t> fast5_idxs;

    std::string read, fast5;
    while(!file.eof()) {
        file >> read >> fast5;
        assert(!read.empty() && !fast5.empty());

        auto itr = fast5_idxs.find(fast5);
        size_t i;
        if (itr == fast5_idxs.end()) {
            i = add_fast5(fast5);
            fast5_idxs[fast5] = i;
        } else {
            i = itr->second;
        }

        add_read(read, i);
    }

    return true;
}

void Fast5Dict::add_read(const std::string &read_id, u32 fast5_idx) {
    reads_[read_id] = fast5_idx;
}

Fast5Read Fast5Dict::operator[](const std::string &read_id) {

    auto itr = reads_.find(read_id);
    if (itr == reads_.end()) {
        return Fast5Read();
    }

    auto i = itr->second;
    if (i != fast5_idx_  && !open_fast5(i)) {
        return Fast5Read();
    }
    fast5_idx_ = i;

    std::string path;
    if (fmt_ == Format::MULTI) {
        path = "/"+MULTI_READ_PREFIX + read_id;
    } else {
        path = get_single_raw_path();
    }

    return Fast5Read(fast5_file_, get_subpaths(path));
}

std::string Fast5Dict::get_read_file(const std::string &read_id) {
    return fast5_files_[reads_[read_id]];
}
