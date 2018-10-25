#ifndef INCL_FMI
#define INCL_FMI

#include <string>
#include "range.hpp"
#include "basepairs.hpp"

class FMI {
    public:

    //virtual void construct(const std::string &seq) = 0;

    virtual void save(const std::string &filename) = 0; 

    virtual Range get_neighbor(Range range, u8 base) const = 0;

    virtual Range get_full_range(u8 base) const = 0;

    virtual u64 sa(u64 i) const = 0;

    virtual u64 size() const = 0;

    bool is_loaded() {return loaded_;}

    protected:
    bool loaded_;
};

#endif
