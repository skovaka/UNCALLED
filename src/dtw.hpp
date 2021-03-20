#ifndef _INCL_DTW
#define _INCL_DTW

#include <vector>
#include <iostream>
#include <cfloat>
#include "util.hpp"
#include "pore_model.hpp"

//TODO don't duplicate klen
//probably should just set constant, or declare elsewhere
const KmerLen DTW_KLEN = KmerLen::k5;

enum class DTWSubSeq {NONE, ROW, COL};
typedef struct {
    DTWSubSeq subseq;
    float dw, hw, vw;
} DTWParams;

const DTWParams 
    DTW_EVENT_GLOB = {
        DTWSubSeq::NONE, 2, 1, 100
    }, DTW_EVENT_QSUB = {
        DTWSubSeq::COL, 2, 1, 100
    }, DTW_EVENT_RSUB = {
        DTWSubSeq::ROW, 2, 1, 100
    }, DTW_RAW_QSUB = {
        DTWSubSeq::COL, 10, 1, 1000
    }, DTW_RAW_RSUB = {
        DTWSubSeq::ROW, 10, 1, 1000
    }, DTW_RAW_GLOB = {
        DTWSubSeq::NONE, 10, 1, 1000
    }, DTW_GLOB = {
        DTWSubSeq::NONE, 1, 1, 1
    };


class DTWp {
    public:

    struct Trace {
        u64 qry, ref;
    };

    DTWp(const std::vector<float> &col_vals,   
        const std::vector<u16> &row_vals,
        const PoreModel<DTW_KLEN> &model,
        const DTWParams &p) : 
        PRMS(p),
        model_(model),
        rvals_(row_vals),
        cvals_(col_vals) {

        mat_.resize(rvals_.size() * cvals_.size());
        bcrumbs_.resize(mat_.size());

        compute_matrix();          
        traceback();
    }

    void compute_matrix() {
        u64 k = 0;
        float cost, ds, hs, vs;
        for (u64 i = 0; i < rvals_.size(); i++) {
            for (u64 j = 0; j < cvals_.size(); j++) {

                cost = costfn(cvals_[j], rvals_[i]);
                ds = dscore(i,j) + (PRMS.dw * cost);
                hs = hscore(i,j) + (PRMS.hw * cost);
                vs = vscore(i,j) + (PRMS.vw * cost);

                if (ds <= hs && ds <= vs) {
                    mat_[k] = ds;
                    bcrumbs_[k++] = Move::D;
                } else if (hs <= vs) {
                    mat_[k] = hs;
                    bcrumbs_[k++] = Move::H;
                } else {
                    mat_[k] = vs;
                    bcrumbs_[k++] = Move::V;
                }
            }
        }
    }

    void traceback() {
        u64 i = rvals_.size()-1, j = cvals_.size()-1;

        switch (PRMS.subseq) {
            case DTWSubSeq::ROW:
                for (u64 k = 0; k < rvals_.size(); k++) {
                    if (mat_[k*cvals_.size() + j] < mat_[i*cvals_.size() + j]) { 
                        i = k;
                    }
                }
                break;
            case DTWSubSeq::COL:
                for (u64 k = 0; k < cvals_.size(); k++) {
                    if (mat_[i*cvals_.size() + k] < mat_[i*cvals_.size() + j]) {
                        j = k;
                    }
                }
                break;
            default:
                break;
        }

        score_sum_ = mat_[i*cvals_.size() + j];

        path_.push_back({j, i});

        u64 k = i*cvals_.size() + j;
        while(!(i == 0 || PRMS.subseq == DTWSubSeq::ROW) ||
              !(j == 0 || PRMS.subseq == DTWSubSeq::COL)) {

            if (i == 0 || bcrumbs_[k] == Move::H) {
                k--;
                j--;
            } else if (j == 0 || bcrumbs_[k] == Move::V) {
                k -= cvals_.size();
                i--;
            } else {
                k -= cvals_.size() + 1;
                i--;
                j--;
            }

            path_.push_back({j, i});
        }
    }

    std::vector<Trace> get_path() {
        return path_;
    }

    float score() {
        return score_sum_;
    }

    float mean_score() {
        return score_sum_ / path_.size();
    }

    protected:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const DTWParams PRMS;
    const PoreModel<DTW_KLEN> model_;

    const std::vector<u16> &rvals_;
    const std::vector<float> &cvals_;

    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector<Trace> path_;
    float score_sum_;

    float costfn(float pA, u16 kmer) {
        return -model_.match_prob(pA,kmer);
    }

    inline float hscore(u64 i, u64 j) {
        if (j > 0) return mat_[cvals_.size()*i + j-1];
        else if (PRMS.subseq == DTWSubSeq::ROW) return 0;
        return MAX_COST;
    }

    inline float vscore(u64 i, u64 j) {
        if (i > 0) return mat_[cvals_.size()*(i-1) + j];
        else if (PRMS.subseq == DTWSubSeq::COL) return 0;
        else return MAX_COST;
    }

    inline float dscore(u64 i, u64 j) {
        if (j > 0 && i > 0) return mat_[cvals_.size()*(i-1) + j-1];
        else if ((j == i) || 
                 (i == 0 && PRMS.subseq == DTWSubSeq::COL) || 
                 (j == 0 && PRMS.subseq == DTWSubSeq::ROW)) return 0;
        return MAX_COST;
    }

    public:

    #ifdef PYBIND
    #define PY_DTW_P_METH(P) c.def(#P, &DTWp::P);
    static void pybind_defs(pybind11::class_<DTWp> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTW_KLEN>,
                             const DTWParams&>());
        PY_DTW_P_METH(get_path)
        PY_DTW_P_METH(score)
        PY_DTW_P_METH(mean_score)
        c.def_property_readonly("path", [](DTWp &d) -> pybind11::array_t<Trace> {
             return pybind11::array_t<Trace>(d.path_.size(), d.path_.data());
        });
        PYBIND11_NUMPY_DTYPE(Trace, qry, ref);
    }
    #endif
};

class DTWd : public DTWp {
    public: 
    DTWd(const std::vector<float> &col_vals,   
        const std::vector<u16> &row_vals,
        const PoreModel<DTW_KLEN> &model,
        const DTWParams &prms) 
        : DTWp(col_vals, row_vals, model, prms) {}

    private:
    float costfn(float pA, u16 kmer) {
        return abs(pA-model_.get_mean(kmer));
    }

    public:
    #ifdef PYBIND
    #define PY_DTW_D_METH(P) c.def(#P, &DTWd::P);
    static void pybind_defs(pybind11::class_<DTWd> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTW_KLEN>,
                             const DTWParams&>());
        PY_DTW_D_METH(get_path)
        PY_DTW_D_METH(score)
        PY_DTW_D_METH(mean_score)
    }
    #endif
};


#endif
