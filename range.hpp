#ifndef INCL_RANGE
#define INCL_RANGE

#include <string>

size_t max(size_t a, size_t b);

size_t min(size_t a, size_t b);

class Range {
    
    public:
    size_t start_, end_; //should be unsigned int, but whatevs

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
