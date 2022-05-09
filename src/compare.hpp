#define COMPARE_HPP

#include "dataframes.hpp"

class Compare {
    public:
    struct Rec {
        int ref;
        float jaccard, mean_ref_dist;
    }

    Compare(AlnCoords a, AlnCoords b) {
        size_t i = 0, j = 0;

        while (i < a.height && j < b.height) {
            

        }

    }

}

#endif
