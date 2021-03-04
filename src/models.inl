#ifndef _INCL_MODELS
#define _INCL_MODELS

#include <unordered_map>
#include "model_r94_dna.inl"
#include "model_r94_rna.inl"

const KmerLen KLEN = KmerLen::k5;

const std::unordered_map<std::string, PoreModel<KLEN>> PORE_MODELS {{
    {"r94_dna_templ",   pmodel_r94_dna_templ},
    {"r94_dna_compl", pmodel_r94_dna_compl},
    {"r94_rna_templ",   pmodel_r94_rna_templ},
    {"r94_rna_compl", pmodel_r94_rna_compl}
}};

#endif
