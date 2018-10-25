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
#define D_MAX_STAY_FRAC 0.5
#define D_MIN_SEED_PROB -3.75
#define D_MIN_MEAN_CONF 6.67
#define D_MIN_TOP_CONF 2

bool load_fast5(const std::string &filename, fast5::File &file) {
    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" << filename << "'\n";
            return false;
        }

        return true;
        
    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        return false;
    }

    return false;
}

void align_fast5s(const std::vector<std::string> &fast5_names,
                 const std::string &bwa_prefix,
                 const std::string &model_fname,
                 const std::string &probfn_fname,
                 u32 seed_len        = D_SEED_LEN,      
                 u32 min_aln_len     = D_MIN_ALN_LEN,
                 u32 min_rep_len     = D_MIN_REP_LEN,
                 u32 max_rep_copy    = D_MAX_REP_COPY,
                 u32 max_consec_stay = D_MAX_CONSEC_STAY,
                 u32 max_paths       = D_MAX_PATHS,
                 float max_stay_frac = D_MAX_STAY_FRAC,
                 float min_seed_prob = D_MIN_SEED_PROB,
                 float min_mean_conf = D_MIN_MEAN_CONF, 
                 float min_top_conf  = D_MIN_TOP_CONF) {

    BwaFMI fmi(bwa_prefix);
    KmerModel model(model_fname, true);

    EventParams event_params = event_detection_defaults;
    event_params.threshold1 = 1.4;
    event_params.threshold2 = 1.1;

    MapperParams aln_params(fmi, model, event_params,
                            probfn_fname, seed_len, min_aln_len,
                            min_rep_len, max_rep_copy, max_consec_stay, 
                            max_paths, max_stay_frac, min_seed_prob,
                            min_mean_conf, min_top_conf);

    Mapper mapper(aln_params);

    for (std::string fast5_name : fast5_names) {
        fast5::File fast5_file;
        if (!load_fast5(fast5_name, fast5_file)) {
            continue;
        }

        auto fast5_info = fast5_file.get_raw_samples_params();
        auto raw_samples = fast5_file.get_raw_samples();
        mapper.new_read(fast5_info.read_id);
        mapper.add_samples(raw_samples);
    }

}

PYBIND11_MODULE(uncalled, m) {
    m.doc() = "UNCALLED";

    m.def("align_fast5s", &align_fast5s, 
          "fast5_names"_a, "bwa_prefix"_a,
          "model_fname"_a, "probfn_fname"_a,
          "seed_len"_a=D_SEED_LEN, "min_aln_len"_a=D_MIN_ALN_LEN,
          "min_rep_len"_a=D_MIN_REP_LEN, "max_rep_copy"_a=D_MAX_REP_COPY,
          "max_consec_stay"_a=D_MAX_CONSEC_STAY, "max_paths"_a=D_MAX_PATHS,
          "max_stay_frac"_a=D_MAX_STAY_FRAC, "min_seed_prob"_a=D_MIN_SEED_PROB,
          "min_mean_conf"_a=D_MIN_MEAN_CONF, "min_top_conf"_a=D_MIN_TOP_CONF);
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
