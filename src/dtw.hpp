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
    float move_cost, stay_cost, skip_cost,
          band_shift;
    i32 band_width;
    std::string band_mode, mm2_paf;
} DtwParams;

const DtwParams 
    DTW_PRMS_DEF = {
        DTWSubSeq::NONE, 1, 1, 1, 0.5, 100, "guided", ""
    }, DTW_PRMS_EVT_GLOB = {
        DTWSubSeq::NONE, 2, 1, 100, 0, 0, "", "",
    };

//   ("method", "guided", str, "DTW method"),
//   ("band_width", 50, int, "DTW band width (only applies to BDTW)"),
//   ("band_shift", 0.5, float, "DTW band shift coefficent (only applies to BDTW)"),
//   ("mm2_paf", None, str, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string."),
//    ("mask_skips", False, bool, "Represent skips as missing data"),

#ifdef PYBIND
#define PY_DTW_PARAM(P, D) p.def_readwrite(#P, &DtwParams::P, D);
void pybind_dtw(py::module_ &m) {
    py::class_<DtwParams> p(m, "DtwParams");
    PY_DTW_PARAM(band_mode, "DTW band mode (\"guided\", \"static\", or \"\"/\"none\")");
    PY_DTW_PARAM(mm2_paf, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string.");
    PY_DTW_PARAM(move_cost, "DTW event move (diagonal) penalty");
    PY_DTW_PARAM(stay_cost, "DTW event stay (horizontal) penalty");
    PY_DTW_PARAM(skip_cost, "DTW event skip (vertical) penalty");
    PY_DTW_PARAM(band_width, "DTW band width");
    PY_DTW_PARAM(band_shift, "DTW band shift");
    //dtwp.def_readwrite("move_cost", &DtwParams::move_cost);
    //dtwp.def_readwrite("stay_cost", &DtwParams::stay_cost);
    //dtwp.def_readwrite("skip_cost", &DtwParams::skip_cost);
    //dtwp.def_readwrite("band_width", &DtwParams::band_width);

    m.attr("DTW_PRMS_DEF") = py::cast(DTW_PRMS_DEF);
    m.attr("DTW_PRMS_EVT_GLOB") = py::cast(DTW_PRMS_EVT_GLOB);
}
#endif

class DTWp {
    public:

    struct Trace {
        u64 qry, ref;
    };

    DTWp(const std::vector<float> &qry_vals,   
        const std::vector<u16> &ref_vals,
        const PoreModel<DTW_KLEN> &model,
        const DtwParams &p) : 
        PRMS(p),
        model_(model),
        ref_vals_(ref_vals),
        qry_vals_(qry_vals) {

        mat_.resize(ref_vals_.size() * qry_vals_.size());
        bcrumbs_.resize(mat_.size());

        compute_matrix();          
        traceback();
    }

    void compute_matrix() {
        u64 k = 0;
        float cost, ds, hs, vs;
        for (u64 r = 0; r < ref_vals_.size(); r++) {
            for (u64 q = 0; q < qry_vals_.size(); q++) {

                cost = costfn(qry_vals_[q], ref_vals_[r]);
                ds = dscore(r,q) + (PRMS.move_cost * cost);
                hs = hscore(r,q) + (PRMS.stay_cost * cost);
                vs = vscore(r,q) + (PRMS.skip_cost * cost);

                if (ds <= hs && ds <= vs) {
                    mat_[k] = ds;
                    bcrumbs_[k] = Move::D;
                } else if (hs <= vs) {
                    mat_[k] = hs;
                    bcrumbs_[k] = Move::H;
                } else {
                    mat_[k] = vs;
                    bcrumbs_[k] = Move::V;
                }

                k++;
            }
        }
    }

    void traceback() {
        u64 r = ref_vals_.size()-1, q = qry_vals_.size()-1;

        switch (PRMS.subseq) {
            case DTWSubSeq::ROW:
                for (u64 k = 0; k < ref_vals_.size(); k++) {
                    if (mat_[k*qry_vals_.size() + q] < mat_[r*qry_vals_.size() + q]) { 
                        r = k;
                    }
                }
                break;
            case DTWSubSeq::COL:
                for (u64 k = 0; k < qry_vals_.size(); k++) {
                    if (mat_[r*qry_vals_.size() + k] < mat_[r*qry_vals_.size() + q]) {
                        q = k;
                    }
                }
                break;
            default:
                break;
        }

        score_sum_ = mat_[r*qry_vals_.size() + q];

        path_.push_back({q, r});

        u64 k = r*qry_vals_.size() + q;
        while(!(r == 0 || PRMS.subseq == DTWSubSeq::ROW) ||
              !(q == 0 || PRMS.subseq == DTWSubSeq::COL)) {

            if (r == 0 || bcrumbs_[k] == Move::H) {
                k--;
                q--;
            } else if (q == 0 || bcrumbs_[k] == Move::V) {
                k -= qry_vals_.size();
                r--;
            } else {
                k -= qry_vals_.size() + 1;
                r--;
                q--;
            }

            path_.push_back({q, r});
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

    //protected:
    public:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const DtwParams PRMS;
    const PoreModel<DTW_KLEN> model_;

    const std::vector<u16> ref_vals_;
    const std::vector<float> qry_vals_;

    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector<Trace> path_;
    float score_sum_;

    float costfn(float pA, u16 kmer) {
        return abs(pA-model_.kmer_current(kmer));
        //return -model_.norm_pdf(pA,kmer);
    }

    inline float hscore(u64 r, u64 q) {
        if (q > 0) return mat_[qry_vals_.size()*r + q-1];
        else if (PRMS.subseq == DTWSubSeq::ROW) return 0;
        return MAX_COST;
    }

    inline float vscore(u64 r, u64 q) {
        if (r > 0) return mat_[qry_vals_.size()*(r-1) + q];
        else if (PRMS.subseq == DTWSubSeq::COL) return 0;
        else return MAX_COST;
    }

    inline float dscore(u64 r, u64 q) {
        if (q > 0 && r > 0) return mat_[qry_vals_.size()*(r-1) + q-1];
        else if ((q == r) || 
                 (r == 0 && PRMS.subseq == DTWSubSeq::COL) || 
                 (q == 0 && PRMS.subseq == DTWSubSeq::ROW)) return 0;
        return MAX_COST;
    }

    public:

    #ifdef PYBIND
    #define PY_DTW_P_METH(P) c.def(#P, &DTWp::P);
    static void pybind_defs(pybind11::class_<DTWp> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTW_KLEN>,
                             const DtwParams&>());
        PY_DTW_P_METH(get_path)
        PY_DTW_P_METH(score)
        PY_DTW_P_METH(mean_score)

        c.def_property_readonly("mat", [](DTWp &d) -> pybind11::array_t<float> {
             return pybind11::array_t<float>(
                     {d.ref_vals_.size(), d.qry_vals_.size()},
                     d.mat_.data()
                     );
        });

        c.def_property_readonly("path", [](DTWp &d) -> pybind11::array_t<Trace> {
             return pybind11::array_t<Trace>(d.path_.size(), d.path_.data());
        });
        PYBIND11_NUMPY_DTYPE(Trace, qry, ref);
    }
    #endif
};

class DTWd : public DTWp {
    public: 
    DTWd(const std::vector<float> &qry_vals,   
        const std::vector<u16> &ref_vals,
        const PoreModel<DTW_KLEN> &model,
        const DtwParams &prms) 
        : DTWp(qry_vals, ref_vals, model, prms) {}

    private:
    float costfn(float pA, u16 kmer) {
        return abs(pA-model_.kmer_current(kmer));
    }

    public:
    #ifdef PYBIND
    #define PY_DTW_D_METH(P) c.def(#P, &DTWd::P);
    static void pybind_defs(pybind11::class_<DTWd, DTWp> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTW_KLEN>,
                             const DtwParams&>());
        PY_DTW_D_METH(get_path)
        PY_DTW_D_METH(score)
        PY_DTW_D_METH(mean_score)
    }
    #endif
};


#endif
