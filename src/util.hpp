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

#ifndef INCL_UTIL
#define INCL_UTIL

#include <string>
#include <fstream>
#include <vector>
#include <cstdint>

#define ALPH_SIZE 4

//Based on github.com/dnbaker/bonsai/blob/master/bonsai/include/util.h
using i8  = std::int8_t;  using u8  = std::uint8_t;
using i16 = std::int16_t; using u16 = std::uint16_t;
using i32 = std::int32_t; using u32 = std::uint32_t;
using i64 = std::int64_t; using u64 = std::uint64_t;

const char BASE_CHARS[] {'A', 'C', 'G', 'T', 'N'};

const u8 BASE_BYTES[] 
     {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //0-15 
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //16-31
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //32-47
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //48-63
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //64-79 (A,C,G)
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //80-95 (T)
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //96-111 (a,c,g)
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};//112-127 (t)

const char BASE_COMP_C[] 
     {'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //0-15  ga
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //16-31 rb
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //32-47 ag
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //48-63 e!
      'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N', //64-79 (A,C,G)
      'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N', //80-95 (T)
      'N','t','N','g','N','N','N','c','N','N','N','N','N','N','N','N', //96-111 (a,c,g)
      'N','N','N','N','a','N','N','N','N','N','N','N','N','N','N','N'};//112-127 (t)

const u8 BASE_COMP_B[] {3, 2, 1, 0};


#endif
