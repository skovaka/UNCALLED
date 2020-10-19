#include <iostream>
#include "event_profiler.hpp"

const EventProfiler::Params EventProfiler::PRMS_DEF {
    win_len        : 25,
    win_stdv_min   : 5,
    win_stdv_range : 3,
    win_mean_range : 4
};
