#ifdef PYBIND
#include "dataframe.hpp"

decltype(AlnCoords::columns) AlnCoords::columns;

#define BIND_DF(name) name::pybind<name>(m, "_"#name);

constexpr AlnCoordsDF::NameArray AlnCoordsDF::names;

void pybind_dataframes(py::module_ &m) {
    AlnCoords::pybind_rec<AlnCoords>(m);

    BIND_DF(AlnCoordsDF)

    PyArray<float>::pybind(m, "PyArrayF32");
    PyArray<int>::pybind(m, "PyArrayI32");  
    //PyArray<u8>::pybind(m, "PyArrayU8");  
    //PyArray<u16>::pybind(m, "PyArrayU16");  
    //PyArray<u32>::pybind(m, "PyArrayU32");  
}

#endif
