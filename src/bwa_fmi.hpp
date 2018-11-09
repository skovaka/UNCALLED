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

#ifndef INCL_BWAFMI
#define INCL_BWAFMI

#include <string>
#include <climits>
#include "util.hpp"
#include "range.hpp"
#include "bwa/bwt.h"
#include "bwa/bntseq.h"

class BwaFMI {
    public:

    BwaFMI();

    BwaFMI(const std::string &prefix);

    Range get_neighbor(Range range, u8 base) const;

    Range get_full_range(u8 base) const;

    u64 sa(u64 i) const;

    u64 size() const;

    u64 translate_loc(u64 sa_loc, std::string &ref_name, u64 &ref_loc) const;

    private:
    bwt_t *index_;
    bntseq_t *bns_;
    bool loaded_;
};

#endif
