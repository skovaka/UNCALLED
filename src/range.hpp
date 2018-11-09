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

#ifndef INCL_RANGE
#define INCL_RANGE

#include <string>

size_t max(size_t a, size_t b);

size_t min(size_t a, size_t b);

class Range {
    
    public:
    //TODO: change to u64
    size_t start_, end_; 

    //Copy constructor
    Range(const Range &prev);

    Range(size_t start, size_t end);

    Range();

    Range& operator=(const Range &r);

    Range split_range(const Range &r);

    Range intersect(const Range &r) const;

    Range merge(const Range &r) const;

    double get_recp_overlap(const Range &r) const;

    bool same_range(const Range &r) const;

    bool intersects(const Range &r) const;

    bool is_valid() const;

    size_t length() const;

    friend bool operator< (const Range &q1, const Range &q2);
};

bool operator< (const Range &q1, const Range &q2);
bool operator== (const Range &q1, const Range &q2);

#endif
