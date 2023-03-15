#include "util.hpp"
    
void pybind_arrays(py::module_ &m) {
    pybind_array<float>(m, "F32");
    pybind_array<bool>(m, "Bool");
    pybind_array<i8> (m, "I8");
    pybind_array<i16>(m, "I16");
    pybind_array<i32>(m, "I32");
    pybind_array<i64>(m, "I64");
    pybind_array<u8> (m, "U8");
    pybind_array<u16>(m, "U16");
    pybind_array<u32>(m, "U32");
    pybind_array<u64>(m, "U64");
}
