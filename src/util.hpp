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

#ifndef _INCL_UTIL
#define _INCL_UTIL

#include <string>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <cassert>

#ifdef PYBIND
#include <pybind11/pybind11.h>
#endif

#ifdef DEBUG_ALL
#define DEBUG_SEEDS
#define DEBUG_PATHS
#define DEBUG_EVENTS
#define DEBUG_THREADS
#define DEBUG_CONFIDENCE
#endif

#if defined(DEBUG_SEEDS) || defined(DEBUG_PATHS) || defined(DEBUG_EVENTS)
#define DEBUG_OUT
#endif

//Based on github.com/dnbaker/bonsai/blob/master/bonsai/include/util.h
using i8  = std::int8_t;  using u8  = std::uint8_t;
using i16 = std::int16_t; using u16 = std::uint16_t;
using i32 = std::int32_t; using u32 = std::uint32_t;
using i64 = std::int64_t; using u64 = std::uint64_t;

class Timer {
    private:
        std::chrono::high_resolution_clock::time_point start;

    public:
        inline Timer() {
            reset();
        }

        inline void reset() {	
            start = std::chrono::high_resolution_clock::now();
        }

        inline double get() {
            return (std::chrono::duration_cast< std::chrono::duration<double> > (std::chrono::high_resolution_clock::now() - start).count()) * 1000.0;
        }

        inline double lap() {
            double ret = get();
            reset();
            return ret;
        }
};

#endif
