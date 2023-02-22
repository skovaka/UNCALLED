#ifdef PYBIND
#include "dataframe.hpp"
#include "util.hpp"

constexpr decltype(AlnCoords::columns) AlnCoords::columns;

#define BIND_DF(name) name::pybind<name>(m, "_"#name);

constexpr AlnCoordsDF::NameArray AlnCoordsDF::names;

void pybind_dataframes(py::module_ &m) {
    AlnCoords::pybind_rec<AlnCoords>(m);

    BIND_DF(AlnCoordsDF)

    pybind_valarray<float>(m, "F32");
    pybind_valarray<i8>(m, "I8");
    pybind_valarray<i16>(m, "I16");
    pybind_valarray<i32>(m, "I32");
    pybind_valarray<i64>(m, "I64");
    pybind_valarray<u8>(m, "U8");
    pybind_valarray<u16>(m, "U16");
    pybind_valarray<u32>(m, "U32");
    pybind_valarray<u64>(m, "U64");

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
