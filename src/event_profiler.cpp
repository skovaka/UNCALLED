#include <iostream>
#include "event_profiler.hpp"

const EventProfiler::Params EventProfiler::PRMS_DEF {
    win_len        : 25,
    win_stdv_min   : 5,
};
