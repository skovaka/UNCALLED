#ifdef PYBIND
#include "dataframe.hpp"

decltype(AlnCoords::columns) AlnCoords::columns;

void pybind_dataframes(py::module_ &m) {
    AlnCoords::pybind_rec<AlnCoords>(m);

    PyArray<float>::pybind(m, "PyArrayF32");
    PyArray<int>::pybind(m, "PyArrayI32");  

}

#endif
