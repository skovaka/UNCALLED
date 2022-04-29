#include "pore_model.hpp"
#include "model_r94_dna.inl"
#include "model_r94_rna.inl"
#include "model_r94_rna_tombo.inl"

//const KmerLen KLEN = 5;

template<>
const PoreModelParams PoreModel<5, u16>::PRMS_DEF {
    name       : "r94_dna",
    reverse    : false,
    complement : false
};

template<>
const PoreModelParams PoreModel<10, u32>::PRMS_DEF {
    name       : "r94_dna",
    reverse    : false,
    complement : false
};

const std::unordered_map<std::string, const std::vector<float> &> 
    PORE_MODEL_PRESETS {{
        {"r94_dna", model_r94_dna_vals},
        {"r94_rna", model_r94_rna_vals},
        {"r94_rna_tombo", model_r94_rna_tombo_vals}
}};

