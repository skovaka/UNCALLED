#include "pore_model.hpp"
#include "model_r94_dna.inl"
#include "model_r94_rna.inl"

const KmerLen KLEN = KmerLen::k5;

template <>
const PoreModel<KLEN>::Params PoreModel<KLEN>::PRMS_DEF {
    name       : "r94_dna",
    reverse    : false,
    complement : false
};

template <>
const std::unordered_map<std::string, const std::vector<float> &> 
    PoreModel<KLEN>::PRESETS {{
        {"r94_dna", model_r94_dna_means_stdvs},
        {"r94_rna", model_r94_rna_means_stdvs}
}};
