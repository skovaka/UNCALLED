#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "fast5.hpp"
#include "mapper.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

#define D_SEED_LEN 22
#define D_MIN_ALN_LEN 25
#define D_MIN_REP_LEN 0
#define D_MAX_REP_COPY 100
#define D_MAX_CONSEC_STAY 8
#define D_MAX_PATHS 10000
#define EVT_WINLEN1 3
#define EVT_WINLEN2 6
#define EVT_THRESH1 1.4
#define EVT_THRESH2 9.0
#define EVT_PEAK_HEIGHT 0.2
#define EVT_MIN_MEAN 30
#define EVT_MAX_MEAN 150
#define D_MAX_STAY_FRAC 0.5
#define D_MIN_SEED_PROB -3.75
#define D_MIN_MEAN_CONF 6.67
#define D_MIN_TOP_CONF 2


PYBIND11_MODULE(uncalled, m) {
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
        .def("align_fast5", 
             &Mapper::align_fast5);
    
}

/*
#include <Python.h>
#include <iostream>
#include "util.hpp"

//sssIIIIIIffff
//
extern "C" {
    static const char *param_names[] = {
        "fast5_list",
        "bwa_prefix",      //s
        "model_fname",     //s
        "probfn_fname",    //s
        "seed_len",        //I
        "min_aln_len",     //I
        "min_repeat_len",  //I
        "max_repeat_copy", //I
        "max_consec_stay", //I
        "max_paths",       //I
        "max_stay_frac",   //f
        "min_seed_prob",   //f
        "min_mean_conf",   //f
        "min_top_conf"     //f
    };
}

static PyObject *uncalled_align_fast5s(PyObject *self,
                                       PyObject *args,
                                       PyObject *keywds) {
    char *fast5_list, *bwa_prefix, *model_fname, *probfn_fname;
    u32 seed_len, min_aln_len, min_rep_len, 
        max_rep_copy, max_consec_stay, max_paths;
    float max_stay_frac, min_seed_prob, min_mean_conf, min_top_conf;
    
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|sssIIIIIIffff", param_names,
            &fast5_list, &bwa_prefix, &model_fname, &probfn_fname,
            &seed_len, &min_aln_len, &min_rep_len, &max_rep_copy, 
            &max_consec_stay, &max_paths, &max_stay_frac, &min_seed_prob, 
            &min_mean_conf, &min_top_conf)) return NULL;

    std::cout << fast5_list << "\n"
              << bwa_prefix << "\n"
              << model_fname << "\n"
              << seed_len << "\n";

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef uncalled_methods[] = {
    {"align_fast5s", (PyCFunction) uncalled_align_fast5s, 
     METH_VARARGS|METH_KEYWORDS, "Align fast5s"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef uncalledmodule = {
    PyModuleDef_HEAD_INIT,
    "uncalled",
    NULL,
    -1,
    uncalled_methods
};

PyMODINIT_FUNC PyInit_uncalled(void) {
    return PyModule_Create(&uncalledmodule);
}
*/
