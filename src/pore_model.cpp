#include "pore_model.hpp"
#include "models/r94_dna.inl"
#include "models/r94_rna.inl"
#include "models/r94_rna_tombo.inl"
#include "models/r9.4_dna_450bps_6mer_npl.inl"

const PoreModelParams PORE_MODEL_PRMS_DEF {
    name       : "r94_dna",
    reverse    : false,
    complement : false,
    k          : 5,
    shift      : 2,
};

using M = PoreModelPreset;
const std::unordered_map<std::string, PoreModelPreset> PORE_MODEL_PRESETS {{
    M{ {"r94_rna", true, false, 5, 2}, model_r94_rna_vals }.map(),
    M{ {"r94_rna_tombo", true, false, 5, 2}, model_r94_rna_tombo_vals }.map(),
    M{ {"r94_dna", false, false, 5, 2}, model_r94_dna_vals }.map(),
    M{ {"r9.4_dna_450bps_6mer_npl", false, false, 6, 2}, model_r9_4_dna_450bps_6mer_npl }.map()
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
