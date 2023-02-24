#ifdef PYBIND
#include "dataframe.hpp"
#include "util.hpp"

constexpr decltype(AlnCoords::columns) AlnCoords::columns;

#define BIND_DF(name) name::pybind<name>(m, "_"#name);

constexpr AlnCoordsDF::NameArray AlnCoordsDF::names;

void pybind_dataframes(py::module_ &m) {
    AlnCoords::pybind_rec<AlnCoords>(m);

    BIND_DF(AlnCoordsDF)
    AlnDF::pybind(m);


    ValArray<float>::pybind(m, "F32");
    ValArray<i8>::pybind(m, "I8");
    ValArray<i16>::pybind(m, "I16");
    ValArray<i32>::pybind(m, "I32");
    ValArray<i64>::pybind(m, "I64");
    ValArray<u8>::pybind(m, "U8");
    ValArray<u16>::pybind(m, "U16");
    ValArray<u32>::pybind(m, "U32");
    ValArray<u64>::pybind(m, "U64");

    //ValArray<float>::pybind(m, "F32");

    IntervalIndex<i64>::pybind(m, "I64");
    IntervalIndex<i32>::pybind(m, "I32");

    PyArray<float>::pybind(m, "PyArrayF32");
    PyArray<int>::pybind(m, "PyArrayI32");  
    PyArray<i64>::pybind(m, "PyArrayI64");  
    PyArray<u8>::pybind(m, "PyArrayU8");  
    PyArray<u16>::pybind(m, "PyArrayU16");  
    PyArray<u32>::pybind(m, "PyArrayU32");  
}

#endif
