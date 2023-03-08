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

#ifndef READ_BUFFER_BC
#define READ_BUFFER_BC

#include <algorithm>
#include "read_buffer.hpp"
#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

class ReadBufferBC : public ReadBuffer {
    public:

    using ReadBuffer::ReadBuffer;
    //ReadBufferBC() : ReadBuffer() {}

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

    void set_moves(const std::vector<u8> &moves, u32 template_start, u32 stride) {
        moves_.assign(moves.begin(), moves.end());
        template_start_ = template_start;
        bce_stride_ = stride;
        bc_loaded_ = true;
    }

    protected:

    bool bc_loaded_ = false;
    u32 template_start_, bce_stride_;
    std::vector<u8> moves_;

    #ifdef PYBIND

    #define PY_READ_BUFFER_BC_METH(N,D) c.def(#N, &ReadBufferBC::N, D);
    #define PY_READ_BUFFER_BC_RPROP(P,D) c.def_property_readonly(#P, &ReadBufferBC::get_##P, D);

    public:

    static void pybind_defs(pybind11::class_<ReadBufferBC, ReadBuffer> &c) {

        c.def(pybind11::init<ReadBufferBC>());
        c.def(pybind11::init<const std::string &, u16, u32, u64, const py::array_t<float> &>());
        c.def("set_moves", &ReadBufferBC::set_moves);
        PY_READ_BUFFER_BC_METH(set_moves, "Set move data");
        PY_READ_BUFFER_BC_RPROP(bc_loaded, "True if basecalling data loaded");
        PY_READ_BUFFER_BC_RPROP(moves, "Guppy BC event moves");
        PY_READ_BUFFER_BC_RPROP(bce_stride, "Guppy BC event length");
        PY_READ_BUFFER_BC_RPROP(template_start, "Sample where guppy basecalling starts");


    }

    #endif
};

#endif
