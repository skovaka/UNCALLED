#include <iostream>
#include "event_profiler.hpp"

const EventProfiler::Params EventProfiler::PRMS_DEF {
    n              : 25,
    win_stdv_min   : 1,
    win_stdv_range : 3,
    win_mean_range : 4
};
