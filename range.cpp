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

    if (end_ > r.end_) {
        start_ = r.end_ + 1;
    } else {
        start_ = 1;
        end_ = 0;
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

double Range::get_recp_overlap(const Range &r) const {
    if (!intersects(r)) {
        return 0;
    }

    return double(intersect(r).length()) / double(merge(r).length());
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
    if (q1.start_ < q2.start_)
        return true;

    if (q1.start_ == q2.start_ && q1.end_ < q2.end_)
        return true;
    
    return false;
}


