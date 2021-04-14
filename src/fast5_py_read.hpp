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

#ifndef PY_FAST5_READ
#define PY_FAST5_READ


#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "read_buffer_bc.hpp"

namespace py = pybind11;

//const std::string
//    GUPPY_SEG_SMRY_PATH = "/Segmentation_000/Summary/segmentation",
//    GUPPY_BC_MOVE_PATH = "/Basecall_1D_000/BaseCalled_template/Move",
//    GUPPY_BC_SMRY_PATH = "/Basecall_1D_000/Summary/basecall_1d_template";

class Fast5PyRead : public ReadBufferBC {
    public:

    Fast5PyRead(py::object read) {
        id_ = read.attr("read_id").cast<std::string>();

        auto handle = read.attr("handle").attr("__getitem__");

        auto global_key = read.attr("global_key").cast<std::string>();
        auto ch_handle = handle(global_key + "channel_id");
        auto ch_attrs = ch_handle.attr("attrs");

        set_channel(
            ch_attrs["channel_number"]
            .attr("astype")("uint16").cast<u16>()
        );

        //Get calibration parameters
        float cal_digit  = ch_attrs["digitisation"].cast<float>(),
              cal_range  = ch_attrs["range"].cast<float>(),
              cal_offset = ch_attrs["offset"].cast<float>(),
              cal_scale = cal_range / cal_digit;

        auto raw_handle = handle(read.attr("raw_dataset_group_name"));
        auto raw_attrs = raw_handle.attr("attrs");
        number_ = raw_attrs["read_number"].cast<u32>();
        start_sample_ = raw_attrs["start_time"].cast<u32>();

        auto raw_data = handle(read.attr("raw_dataset_name"));


        u32 raw_len = raw_data.attr("size").cast<u32>(); //TODO use request() info?

        //Write into ReadBuffer method, don't repeat
        u32 min_sample = std::min(PRMS.start_chunk * PRMS.chunk_len(), raw_len);
        if (PRMS.skip_notempl && min_sample < template_start_) {
            min_sample = template_start_;
        }
        start_sample_ += min_sample; //Update global sample start
        u32 max_sample = std::min(min_sample + (PRMS.max_chunks * PRMS.chunk_len()), raw_len);
        full_duration_ = max_sample-min_sample;//TODO full_duration_ is confusing
        chunk_count_ = (full_duration_ / PRMS.chunk_len()) + (full_duration_ % PRMS.chunk_len() != 0);

        signal_.reserve(max_sample-min_sample);
        
        //TODO error checks based on https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#vectorizing-functions ?
        //TODO support other formats than u16? like live chunk loading or whatever?

        auto raw_slice = raw_data.attr("__getitem__")
                         (py::slice(min_sample, max_sample, 1));

        //auto raw_buf = py::buffer(raw_data.attr("__getitem__")(py::tuple()));
        auto raw_buf = py::buffer(raw_slice);
        auto raw_info = raw_buf.request();
        auto raw_arr = static_cast<u16 *>(raw_info.ptr);

        //for (size_t i = min_sample; i < max_sample; i++) {
        for (size_t i = 0; i < raw_info.shape[0]; i++) {
            //TODO this is wrong? add offset before scaling
            signal_.push_back(cal_scale * (raw_arr[i] + cal_offset));
        }

        //py::print(id_, get_channel(), raw_buf.attr("__class__"));
    }

    private:

    #ifdef PYBIND

    public:

    static void pybind_defs(pybind11::class_<Fast5PyRead, ReadBufferBC> &c) {
        c.def(pybind11::init<py::object>());
    }

    #endif
};

#endif
