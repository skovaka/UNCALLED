#ifdef PYBIND
#include "dataframe.hpp"

#define BIND_DF(name) name::pybind<name>(m, #name);

constexpr AlnCoords::NameArray AlnCoords::names;

void pybind_dataframes(py::module_ &m) {
    PyArray<float>::pybind(m, "F32");
    PyArray<int>::pybind(m, "I32");  

    BIND_DF(AlnCoords);
}

#endif
