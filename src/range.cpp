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

#include "range.hpp"

size_t max(size_t a, size_t b) {
    return a > b ? a : b;
}

size_t min(size_t a, size_t b) {
    return a < b ? a : b;
}

Range::Range(const Range &prev)
    : start_(prev.start_), 
      end_(prev.end_) {}

Range::Range(size_t start, size_t end) : start_(start), end_(end) {}


Range::Range() : start_(1), end_(0) {}

bool Range::intersects(const Range &q) const {
    return is_valid() && q.is_valid() &&
           !(start_ > q.end_ || end_ < q.start_) &&
           !(q.start_ > end_ || q.end_ < start_);
}

size_t Range::length() const {
    return end_ - start_ + 1;
}

Range Range::split_range(const Range &r) { 

    Range left;
    if (start_ < r.start_) {
        left = Range(*this);
        left.end_ = r.start_ - 1;
    }

    if (start_ <= r.end_) {
        if (end_ > r.end_) {
            start_ = r.end_ + 1;
        } else {
            start_ = 1;
            end_ = 0;
        }
    }

    return left;
}

Range Range::intersect(const Range &r) const {
    if (!intersects(r)) {
        return Range();
    }
    
    return Range(max(start_, r.start_), min(end_, r.end_));
}

Range Range::merge(const Range &r) const {
    if (!intersects(r)) {
        return Range();
    }

    return Range(min(start_, r.start_), max(end_, r.end_));
}

float Range::get_recp_overlap(const Range &r) const {
    if (!intersects(r)) {
        return 0;
    }

    return float(intersect(r).length()) / float(merge(r).length());
}

Range& Range::operator=(const Range& prev) {
    start_ = prev.start_;
    end_ = prev.end_;
    return *this;
}


bool Range::same_range(const Range &q) const {
    return start_ == q.start_ && end_ == q.end_;
}

bool Range::is_valid() const {
    return start_ <= end_;
}


bool operator< (const Range &q1, const Range &q2) {
    return q1.start_ < q2.start_ ||
           (q1.start_ == q2.start_ && q1.end_ < q2.end_);
}


bool operator== (const Range &q1, const Range &q2) {
    return q1.start_ == q2.start_ && q1.end_ == q2.end_;
}
