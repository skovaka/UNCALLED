#ifdef PYBIND
#include "dataframe.hpp"

#define BIND_DF(name) name::pybind<name>(m, "_"#name);

constexpr AlnCoordsOld::NameArray AlnCoordsOld::names;

void pybind_dataframes(py::module_ &m) {
    PYBIND11_NUMPY_DTYPE(AlnCoord, ref, start, end);
    PyArray<AlnCoord>::pybind(m, "AlnCoordsNew");

    PyArray<float>::pybind(m, "PyArrayF32");
    PyArray<int>::pybind(m, "PyArrayI32");  

    BIND_DF(AlnCoordsOld);
}

#endif
