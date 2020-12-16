#ifndef _INCL_DTW
#define _INCL_DTW

#include <vector>
#include <iostream>
#include <cfloat>
#include "util.hpp"

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
    };


template < class ColT, class RowT, typename Func >
class DTW {
    public:

    DTW(const std::vector<ColT> &col_vals,   
        const std::vector<RowT> &row_vals,
        const DTWParams &p,
        const Func &fn) : 
        PRMS(p),
        fn_(fn),
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

                cost = fn_(rvals_[i], cvals_[j]);
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

        path_.push_back(std::pair<u64, u64>(j, i));

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

            path_.push_back(std::pair<u64, u64>(j, i));
        }
    }

    std::vector< std::pair<u64, u64> > get_path() {
        return path_;
    }

    float score() {
        return score_sum_;
    }

    float mean_score() {
        return score_sum_ / path_.size();
    }

    //friend std::ostream &operator<< (std::ostream &out, const DTW &a);
    void print_path(std::ostream &out) {
        for (auto p = path_.rbegin(); p != path_.rend(); p++) {
            out << p->first << "\t"
                << p->second << "\t"
                << rvals_[p->first] << "\t"
                << cvals_[p->second] << "\t"
                << fn_(rvals_[p->first], cvals_[p->second]) << "\t"
                << "\n";
        }
    } 

    private:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const DTWParams PRMS;
    const Func &fn_;

    const std::vector<RowT> &rvals_;
    const std::vector<ColT> &cvals_;

    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector< std::pair<u64, u64> > path_;
    float score_sum_;

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

    #ifdef PYBIND
    
    static void pybind_defs(pybind11::class_<DTW> &c) {
    }
    #endif
};

float dtwcost_r94p(u16 k, float e) {
    return -pmodel_r94_template.match_prob(e,k);
}
class DTWr94p : public DTW<float, u16, decltype(dtwcost_r94p)> {
    public:
    DTWr94p(const std::vector<float> &means,
            const std::vector<u16> &kmers,
            const DTWParams &prms) 
        : DTW(means, kmers, prms, dtwcost_r94p) {}

    #ifdef PYBIND
    #define PY_DTW_R94P_METH(P) c.def(#P, &DTWr94p::P);
    static void pybind_defs(pybind11::class_<DTWr94p> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const DTWParams&>());
        PY_DTW_R94P_METH(get_path)
        PY_DTW_R94P_METH(score)
        PY_DTW_R94P_METH(mean_score)
    }
    #endif
    
};

float dtwcost_r94d(u16 k, float e) {
    return abs(e-pmodel_r94_template.get_mean(k));
}
class DTWr94d : public DTW<float, u16, decltype(dtwcost_r94d)> {
    public:
    DTWr94d(const std::vector<float> &means,
            const std::vector<u16> &kmers,
            const DTWParams &prms) 
        : DTW(means, kmers, prms, dtwcost_r94d) {}

    #ifdef PYBIND
    #define PY_DTW_R94D_METH(P) c.def(#P, &DTWr94d::P);
    static void pybind_defs(pybind11::class_<DTWr94d> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const DTWParams&>());
        PY_DTW_R94D_METH(get_path)
        PY_DTW_R94D_METH(score)
        PY_DTW_R94D_METH(mean_score)
    }
    #endif
};

#endif
