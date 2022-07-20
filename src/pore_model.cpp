#include "pore_model.hpp"
#include "models/r94_dna.inl"
#include "models/r94_rna_new.inl"
#include "models/r94_rna_tombo_new.inl"

const PoreModelParams PORE_MODEL_PRMS_DEF {
    name       : "r94_dna",
    reverse    : false,
    complement : false,
    k          : 5,
    shift      : 2,
};

const std::unordered_map<std::string, PoreModelPreset> PORE_MODEL_PRESETS {{
    model_r94_dna.map(), model_r94_rna.map(), model_r94_rna_tombo.map(),
}};

#ifdef PYBIND

void pybind_pore_model_params(py::module_ &m) {
    py::class_<PoreModelParams> p(m, "PoreModelParams");
    p.def(py::init<>());
    p.def(py::init<PoreModelParams>());
    PY_MODEL_PARAM(name, "Model preset name or TSV filename");
    PY_MODEL_PARAM(k, "K-mer length");
    PY_MODEL_PARAM(shift, "K-mer shift");
    PY_MODEL_PARAM(reverse, "Will reverse (flip) k-mer sequences if True");
    PY_MODEL_PARAM(complement, "Will complement k-mer sequences if True");

    py::class_<PoreModelPreset> pre(m, "PoreModelPresets");
    pre.def_readwrite("prms", &PoreModelPreset::prms);
    //pre.def_readwrite("vals", &PoreModelPreset::vals);

    m.attr("PORE_MODEL_PRESETS") = py::cast(PORE_MODEL_PRESETS);
}
#endif
