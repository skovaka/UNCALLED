#ifndef INCL_SDSLFMI
#define INCL_SDSLFMI

#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <climits>
#include "fmi.hpp"
#include "basepairs.hpp"
#include "range.hpp"

class SdslFMI : public FMI {
    public:

    SdslFMI();

    SdslFMI(const std::string &filename);

    void construct(const std::string &seq);

    void save(const std::string &filename); 

    Range get_neighbor(Range range, base_t base) const;

    Range get_full_range(base_t base) const;

    size_t sa(size_t i) const;

    size_t size() const;

    private:
    sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector>, 64, UINT_MAX> index_;
};

#endif
