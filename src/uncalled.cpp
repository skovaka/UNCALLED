#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "self_align_ref.hpp"
#include "mapper.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(align, m) {
    m.doc() = "UNCALLED";

    py::class_<MapperParams>(m, "MapperParams")
        .def(py::init<const std::string &,
                      const std::string &,
                      const std::string &,
                      u32, u32, u32, u32, u32, u32,  u32, u32, u32,
                      float, float, float, float, float, 
                      float, float, float, float>());

    py::class_<Mapper>(m, "Mapper")
        .def(py::init<MapperParams &>())
        .def("map_fast5", 
             &Mapper::map_fast5);
    
    m.def("self_align", &self_align);
}

