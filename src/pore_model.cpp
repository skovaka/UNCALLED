#include "pore_model.hpp"
#include "models/r94_dna.inl"
#include "models/r94_rna.inl"
#include "models/r94_rna_tombo.inl"
#include "models/r9.4_dna_450bps_6mer_npl.inl"

const PoreModelParams PORE_MODEL_PRMS_DEF {
    name          : "",
    k             : -1,
    shift         : -1,
    pa_mean       : -1,
    pa_stdv       : -1,
    norm_max      : 5.0,
    sample_rate   : 4000, 
    bases_per_sec : 450,
    reverse       : false,
    complement    : false,
    flowcell      : "",
    kit           : "",
};

#ifdef PYBIND

void pybind_pore_model_params(py::module_ &m) {
    py::class_<PoreModelParams> p(m, "PoreModelParams");
    p.def(py::init<>());
    p.def(py::init<PoreModelParams>());
    p.def("to_key", 
        [&](const PoreModelParams &p) 
            -> std::tuple<std::string, bool, bool> {
                return {p.name, p.reverse, p.complement};
            });
    p.def("norm_to_inorm", py::vectorize(static_cast<inorm_t (PoreModelParams::*) (float) const>(&PoreModelParams::norm_to_inorm)));
    p.def("inorm_to_norm", py::vectorize(static_cast<float (PoreModelParams::*) (inorm_t) const>(&PoreModelParams::inorm_to_norm)));
    p.def("norm_to_pa", py::vectorize(static_cast<float (PoreModelParams::*) (float) const>(&PoreModelParams::norm_to_pa)));
    p.def("pa_to_norm", py::vectorize(static_cast<float (PoreModelParams::*) (float) const>(&PoreModelParams::pa_to_norm)));
    p.def_property_readonly("inorm_scale", &PoreModelParams::inorm_scale);
    PY_MODEL_PARAM(name, "Model preset name or TSV filename");
    PY_MODEL_PARAM(k, "K-mer length");
    PY_MODEL_PARAM(shift, "K-mer shift");
    PY_MODEL_PARAM(norm_max, "K-mer shift");
    PY_MODEL_PARAM(sample_rate, "Raw signal samples per second");
    PY_MODEL_PARAM(bases_per_sec, "Average bases sequenced per second (approximate)");
    PY_MODEL_PARAM(pa_mean, "Current mean picoamp mean");
    PY_MODEL_PARAM(pa_stdv, "Current mean picoamp stdv");
    PY_MODEL_PARAM(reverse, "Will reverse (flip) k-mer sequences if True");
    PY_MODEL_PARAM(complement, "Will complement k-mer sequences if True");
    PY_MODEL_PARAM(flowcell, "Flowcell used for sequencing (e.g. FLO-MIN106)");
    PY_MODEL_PARAM(kit, "Kit used for sequencing (e.g. SQK-LSK109)");
    p.def("has_workflow", &PoreModelParams::has_workflow);
    p.def("get_workflow", &PoreModelParams::get_workflow);

    //py::class_<PoreModelPreset> pre(m, "PoreModelPresets");
    //pre.def_readwrite("prms", &PoreModelPreset::prms);
    //pre.def_readwrite("vals", &PoreModelPreset::vals);

    //m.attr("PORE_MODEL_PRESETS") = py::cast(PORE_MODEL_PRESETS);
}
#endif
