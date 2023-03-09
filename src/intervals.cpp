#ifdef PYBIND
#include "intervals.hpp"
#include "util.hpp"

void pybind_intervals(py::module_ &m) {
    IntervalIndex<i64>::pybind(m, "I64");
    IntervalIndex<i32>::pybind(m, "I32");
}

#endif
