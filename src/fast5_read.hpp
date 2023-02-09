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

#include <algorithm>
#include "hdf5_tools.hpp"
#include "read_buffer_bc.hpp"

const std::string GUPPY_SEG_SMRY = "/Summary/segmentation",
                  GUPPY_BC_MOVE = "/BaseCalled_template/Move",
                  GUPPY_BC_SMRY = "/Summary/basecall_1d_template";

class Fast5Read : public ReadBufferBC {
    public:

    struct Paths {
        std::string raw, channel, basecall, segmentation;
    };

    Fast5Read() : ReadBufferBC() {}

    Fast5Read(const hdf5_tools::File &file, 
              const Paths &paths) {
              //const std::string &raw_path, 
              //const std::string &ch_path, 
              //const std::string &anl_path) {

        filename_ = file.file_name();

        //Get read ID (name), number, and start_time
        auto raw_attrs = file.get_attr_map(paths.raw);
        id_ = raw_attrs["read_id"];
        number_ = atoi(raw_attrs["read_number"].c_str());
        start_sample_ = atoi(raw_attrs["start_time"].c_str());

        //Get read channel
        auto ch_attrs = file.get_attr_map(paths.channel);
        channel_idx_ = atoi(ch_attrs["channel_number"].c_str()) - 1;

        //And callibration data
        float cal_digit = atof(ch_attrs["digitisation"].c_str()),
              cal_range = atof(ch_attrs["range"].c_str()),
              cal_offset = atof(ch_attrs["offset"].c_str());

        //Load BC data if path exists
        if (!paths.basecall.empty() && file.group_exists(paths.basecall)) {
            bc_loaded_ = true;

            //Template start from guppy segmentation
            std::string seg_path = paths.segmentation + GUPPY_SEG_SMRY;
            auto seg_attrs = file.get_attr_map(seg_path);
            template_start_ = atoi(seg_attrs["first_sample_template"].c_str());
    
            //Guppy BaseCalled Event constant length (stride)
            auto bc_attrs = file.get_attr_map(paths.basecall + GUPPY_BC_SMRY);
            bce_stride_ = atoi(bc_attrs["block_stride"].c_str());

            //Read guppy bce moves
            std::string bc_move_path = paths.basecall + GUPPY_BC_MOVE;
            file.read(bc_move_path, moves_);

        } else {
            bc_loaded_ = false;
            template_start_ = 0;
        }

        std::string sig_path = paths.raw + "/Signal";
        std::vector<i16> int_data; 
        file.read(sig_path, int_data);

        u32 raw_len = static_cast<u32>(int_data.size());

        //Determine min sample
        //
        u32 min_sample = std::min(PRMS.start_chunk * PRMS.chunk_len(), raw_len);
        if (PRMS.skip_notempl && min_sample < template_start_) {
            min_sample = template_start_;
        }
        start_sample_ += min_sample; //Update global sample start

        //Set read length
        u32 max_sample = std::min(min_sample + (PRMS.max_chunks * PRMS.chunk_len()), raw_len);
        full_duration_ = max_sample-min_sample;
        chunk_count_ = (full_duration_ / PRMS.chunk_len()) + (full_duration_ % PRMS.chunk_len() != 0);

        //Fill in calibrated signal
        signal_.reserve(full_duration_);
        for (u32 i = min_sample; i < max_sample; i++) {
            float calibrated = (cal_range * (int_data[i] + cal_offset) / cal_digit);
            signal_.push_back(calibrated);
        }
        assert(full_duration_ == signal_.size());

        //Trim BC moves

        if (bc_loaded_) {
            u32 min_bce = 0;
            if (min_sample > template_start_) {
                min_bce = (min_sample - template_start_) / bce_stride_;
            }

            u32 max_bce = (max_sample - template_start_ + 1) / bce_stride_;

            if (max_bce - min_bce < moves_.size()) {
                moves_ = std::vector<u8>(&moves_[min_bce], &moves_[max_bce-1]);
            }
        }

    }

    protected:
    //std::string filename;

    #ifdef PYBIND

    #define PY_FAST5_ATTR(P,D) c.def_readonly(#P, &Fast5Read::P, D);

    public:

    static void pybind_defs(pybind11::class_<Fast5Read, ReadBufferBC> &c) {

        c.def("__repr__", [](Fast5Read &r) -> std::string {
            return "<Fast5Read id=" + r.get_id() + ">"; //TODO add channel, number, length
        });


    }

    #endif
};
