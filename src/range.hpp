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

#ifndef _INCL_RANGE
#define _INCL_RANGE

#include <string>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif 

i64 max(i64 a, i64 b);

i64 min(i64 a, i64 b);


class Range {
    
    public:
    i64 start_, end_; 

    //Copy constructor
    Range(const Range &prev);

    Range(i64 start, i64 end);

    Range();

    Range& operator=(const Range &r);

    Range split_range(const Range &r);

    Range intersect(const Range &r) const;

    Range merge(const Range &r) const;

    float get_recp_overlap(const Range &r) const;

    bool intersects(const Range &r) const;

    bool is_valid() const;

    i64 length() const;

    friend bool operator< (const Range &q1, const Range &q2);
    friend bool operator== (const Range &q1, const Range &q2);

    #ifdef PYBIND

    #define PY_RANGE_METH(P) c.def(#P, &Range::P);

    static void pybind_defs(pybind11::class_<Range> &c) {
        c.def(pybind11::init<i64, i64>());
        PY_RANGE_METH(intersect)
        PY_RANGE_METH(merge)
        PY_RANGE_METH(get_recp_overlap)
        PY_RANGE_METH(intersects)
        PY_RANGE_METH(is_valid)
        PY_RANGE_METH(length)
        c.def_property_readonly("is_valid", &Range::is_valid);
        c.def_property_readonly("length", &Range::length);
        c.def_readwrite("start", &Range::start_);
        c.def_readwrite("end", &Range::end_);

        //PYBIND11_NUMPY_DTYPE(Range, start, end);
    }
    #endif

};

bool operator< (const Range &q1, const Range &q2);
bool operator== (const Range &q1, const Range &q2);

#endif
