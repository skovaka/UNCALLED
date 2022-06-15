#include "pore_model.hpp"
#include "model_r94_dna.inl"
#include "model_r94_rna.inl"
#include "model_r94_rna_tombo.inl"

//const KmerLen KLEN = 5;

//template<>
//const PoreModelParams PoreModel<5, u16>::PRMS_DEF {
//    name       : "r94_dna",
//    reverse    : false,
//    complement : false
//};
//
//template<>
//const PoreModelParams PoreModel<10, u32>::PRMS_DEF {
//    name       : "r94_dna",
//    reverse    : false,
//    complement : false
//};

const PoreModelParams PORE_MODEL_PRMS_DEF {
    name       : "r94_dna",
    reverse    : false,
    complement : false
};

template<>
const PresetMap 
    PoreModel<5>::PRESETS {{
            {"r94_dna", model_r94_dna_vals},
            {"r94_rna", model_r94_rna_vals},
            {"r94_rna_tombo", model_r94_rna_tombo_vals}
}}; 
 
template<> const PresetMap PoreModel<4>::PRESETS = {};
template<> const PresetMap PoreModel<6>::PRESETS = {};
template<> const PresetMap PoreModel<7>::PRESETS = {};
template<> const PresetMap PoreModel<8>::PRESETS = {};
template<> const PresetMap PoreModel<9>::PRESETS = {};
template<> const PresetMap PoreModel<10>::PRESETS = {};
template<> const PresetMap PoreModel<11>::PRESETS = {};
template<> const PresetMap PoreModel<12>::PRESETS = {};

template <> const KmerLen PoreModel<4>::SHIFT = 1;
template <> const KmerLen PoreModel<5>::SHIFT = 2;
template <> const KmerLen PoreModel<6>::SHIFT = 2;
template <> const KmerLen PoreModel<7>::SHIFT = 3;
template <> const KmerLen PoreModel<8>::SHIFT = 3;
template <> const KmerLen PoreModel<9>::SHIFT = 4;
template <> const KmerLen PoreModel<10>::SHIFT = 4;
template <> const KmerLen PoreModel<12>::SHIFT = 5;

#ifdef PYBIND

void pybind_pore_model_params(py::module_ &m) {
    py::class_<PoreModelParams> p(m, "PoreModelParams");
    p.def(py::init<>());
    p.def(py::init<PoreModelParams>());
    PY_MODEL_PARAM(name, "Model preset name or TSV filename");
    PY_MODEL_PARAM(reverse, "Will reverse (flip) k-mer sequences if True");
    PY_MODEL_PARAM(complement, "Will complement k-mer sequences if True");
}
#endif
