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

#include "hdf5_tools.hpp"
#include "read_buffer.hpp"

const std::string
    GUPPY_SEG_SMRY_PATH = "/Segmentation_000/Summary/segmentation",
    GUPPY_BC_MOVE_PATH = "/Basecall_1D_000/BaseCalled_template/Move",
    GUPPY_BC_SMRY_PATH = "/Basecall_1D_000/Summary/basecall_1d_template";

class Fast5Read : public ReadBuffer {
    public:

    Fast5Read(const hdf5_tools::File &file, 
              const std::string &raw_path, 
              const std::string &ch_path, 
              const std::string &anl_path) {

        for (auto a : file.get_attr_map(raw_path)) {
            if (a.first == "read_id") {
                id_ = a.second;
            } else if (a.first == "read_number") {
                number_ = atoi(a.second.c_str());
            } else if (a.first == "start_time") {
                start_sample_ = atoi(a.second.c_str());
            }
        }

        float cal_digit = 1, cal_range = 1, cal_offset = 0;
        for (auto a : file.get_attr_map(ch_path)) {
            if (a.first == "channel_number") {
                channel_idx_ = atoi(a.second.c_str()) - 1;
            } else if (a.first == "digitisation") {
                cal_digit = atof(a.second.c_str());
            } else if (a.first == "range") {
                cal_range = atof(a.second.c_str());
            } else if (a.first == "offset") {
                cal_offset = atof(a.second.c_str());
            }
        }

        template_start_ = 0;

        if (!anl_path.empty() && file.group_exists(anl_path)) {
            bc_loaded_ = true;

            //Template start from guppy segmentation
            std::string seg_path = anl_path + GUPPY_SEG_SMRY_PATH;
            for (auto a : file.get_attr_map(seg_path)) {
                if (a.first == "first_sample_template") {
                    template_start_ = atoi(a.second.c_str());
                }
            }
    
            //Guppy BaseCalled Event constant length (stride)
            auto bc_attrs = file.get_attr_map(anl_path + GUPPY_BC_SMRY_PATH);
            bce_stride_ = atoi(bc_attrs["block_stride"].c_str());

            //Read guppy bce moves
            std::string bc_move_path = anl_path + GUPPY_BC_MOVE_PATH;
            file.read(bc_move_path, moves_);

        } else {
            bc_loaded_ = false;
        }

        //TODO clean up constructors
        //don't use inheritance

        std::string sig_path = raw_path + "/Signal";
        std::vector<i16> int_data; 
        file.read(sig_path, int_data);

        chunk_count_ = (int_data.size() / PRMS.chunk_len()) + (int_data.size() % PRMS.chunk_len() != 0);

        u32 start_sample = PRMS.start_chunk * PRMS.chunk_len();

        if (PRMS.skip_notempl && start_sample < template_start_) {
            start_sample = template_start_;
        }

        if (int_data.size() <= start_sample) {
            int_data.clear();
        } else {
            int_data = std::vector<i16>(&int_data[start_sample], &int_data[int_data.size()]);
        }

        if (chunk_count_ > PRMS.max_chunks) {
            chunk_count_ = PRMS.max_chunks;
            int_data.resize(chunk_count_ * PRMS.chunk_len());
        }

        //signal_.reserve(int_data.size());

        //signal_.assign(int_data.begin(), int_data.end());

        //for (u64 i = template_start_; i < int_data.size(); i++) {
        for (u16 raw : int_data) {
            float calibrated = (cal_range * raw / cal_digit) + cal_offset;
            signal_.push_back(calibrated);
        }

        full_duration_ = signal_.size();
    }

    u32 get_template_start() const {
        return template_start_;
    }

    const std::vector<u8> &get_moves() const {
        return moves_;
    }

    bool get_bc_loaded() const {
        return bc_loaded_;
    }

    u32 get_bce_stride() const {
        return bce_stride_;
    }

    private:

    bool bc_loaded_;
    u32 template_start_, bce_stride_;
    std::vector<u8> moves_;

    #ifdef PYBIND

    #define PY_FAST5_READ_METH(N,D) c.def(#N, &Fast5Read::N, D);
    #define PY_FAST5_READ_RPROP(P,D) c.def_property_readonly(#P, &Fast5Read::get_##P, D);

    public:

    static void pybind_defs(pybind11::class_<Fast5Read, ReadBuffer> &c) {

        //c.def(pybind11::init<ReadBuffer>());
        PY_FAST5_READ_RPROP(bc_loaded, "True if basecalling data loaded");
        PY_FAST5_READ_RPROP(moves, "Guppy BC event moves");
        PY_FAST5_READ_RPROP(bce_stride, "Guppy BC event length");
        PY_FAST5_READ_RPROP(template_start, "Sample where guppy basecalling starts");


    }

    #endif
};
