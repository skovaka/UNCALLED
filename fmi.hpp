#ifndef INCL_FMI
#define INCL_FMI

#include <string>
#include "range.hpp"
#include "basepairs.hpp"

class FMI {
    public:

    //virtual void construct(const std::string &seq) = 0;

    virtual void save(const std::string &filename) = 0; 

    virtual Range get_neighbor(Range range, base_t base) const = 0;

    virtual Range get_full_range(base_t base) const = 0;

    virtual size_t sa(size_t i) const = 0;

    virtual size_t size() const = 0;

    bool is_loaded() {return loaded_;}
    protected:
    //private:
    bool loaded_;
};

#endif
