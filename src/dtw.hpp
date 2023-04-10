#ifndef _INCL_DTW
#define _INCL_DTW

#include <vector>
#include <iostream>
#include <cfloat>
#include "util.hpp"
#include "pore_model.hpp"
#include "aln.hpp"

enum class DTWSubSeq {NONE, ROW, COL};
enum class DTWCostFn {ABS_DIFF, NORM_PDF, Z_SCORE};

struct DtwParams {
    DTWSubSeq subseq;
    float move_cost, stay_cost, skip_cost,
          band_shift;
    i32 del_max, ins_max, band_width, iterations;
    std::string norm_mode, band_mode, cost_fn;
    bool save_bands;
};

extern const DtwParams DTW_PRMS_DEF, DTW_PRMS_EVT_GLOB;

#ifdef PYBIND
#define PY_DTW_PARAM(P, D) p.def_readwrite(#P, &DtwParams::P, D);
void pybind_dtw(py::module_ &m);
#endif


template <typename ModelType>
class GlobalDTW {
    public:

    using KmerType = typename ModelType::kmer_t;

    struct Trace {
        u64 qry, ref;
    };

    GlobalDTW(const std::vector<float> &qry_vals,   
        const std::vector<KmerType> &ref_vals,
        const ModelType &model,
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
            //for (u64 q = event_start_; q < event_end_; q++) {

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

    size_t ref_size() const {
        return (ref_vals_.size());
    }

    size_t qry_size() const {
        return (qry_vals_.size());
    }

    void traceback() {
        u64 r = ref_vals_.size()-1, q = qry_size()-1;//qry_vals_.size()-1;

        switch (PRMS.subseq) {
            case DTWSubSeq::ROW:
                for (u64 k = 0; k < ref_vals_.size(); k++) {
                    if (mat_[k*qry_size() + q] < mat_[r*qry_size() + q]) { 
                        r = k;
                    }
                }
                break;
            case DTWSubSeq::COL:
                for (u64 k = 0; k < qry_size(); k++) {
                    if (mat_[r*qry_size() + k] < mat_[r*qry_size() + q]) {
                        q = k;
                    }
                }
                break;
            default:
                break;
        }

        score_sum_ = mat_[r*qry_size() + q];

        path_.push_back({q, r});

        u64 k = r*qry_size() + q;
        while(!(r == 0 || PRMS.subseq == DTWSubSeq::ROW) ||
              !(q == 0 || PRMS.subseq == DTWSubSeq::COL)) {

            if (r == 0 || bcrumbs_[k] == Move::H) {
                k--;
                q--;
            } else if (q == 0 || bcrumbs_[k] == Move::V) {
                k -= qry_size();
                r--;
            } else {
                k -= qry_size() + 1;
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
    const ModelType model_;

    const std::vector<KmerType> ref_vals_;
    const std::vector<float> qry_vals_;

    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector<Trace> path_;
    float score_sum_;

    float costfn(float pA, KmerType kmer) {
        return abs(pA-model_.kmer_current(kmer));
        //return -model_.norm_pdf(pA,kmer);
    }

    inline float hscore(u64 r, u64 q) {
        if (q > 0) return mat_[qry_size()*r + q-1];
        else if (PRMS.subseq == DTWSubSeq::ROW) return 0;
        return MAX_COST;
    }

    inline float vscore(u64 r, u64 q) {
        if (r > 0) return mat_[qry_size()*(r-1) + q];
        else if (PRMS.subseq == DTWSubSeq::COL) return 0;
        else return MAX_COST;
    }

    inline float dscore(u64 r, u64 q) {
        if (q > 0 && r > 0) return mat_[qry_size()*(r-1) + q-1];
        else if ((q == r) || 
                 (r == 0 && PRMS.subseq == DTWSubSeq::COL) || 
                 (q == 0 && PRMS.subseq == DTWSubSeq::ROW)) return 0;
        return MAX_COST;
    }

    public:

    #ifdef PYBIND
    #define PY_DTW_P_METH(P) c.def(#P, &GlobalDTW<ModelType>::P);
    static void pybind_defs(pybind11::module_ &m, const std::string &suffix) {
        pybind11::class_<GlobalDTW> c(m, ("GlobalDTW" + suffix).c_str());
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<KmerType>&,
                             const ModelType&,
                             const DtwParams&>());
        PY_DTW_P_METH(get_path)
        PY_DTW_P_METH(score)
        PY_DTW_P_METH(mean_score)

        c.def_property_readonly("mat", [](GlobalDTW<ModelType> &d) -> pybind11::array_t<float> {
             return pybind11::array_t<float>(
                     {d.ref_vals_.size(), d.qry_size()},
                     d.mat_.data()
                     );
        });

        c.def_property_readonly("path", [](GlobalDTW<ModelType> &d) -> pybind11::array_t<Trace> {
             return pybind11::array_t<Trace>(d.path_.size(), d.path_.data());
        });
        //PYBIND11_NUMPY_DTYPE(Trace, qry, ref);
    }
    #endif
};

struct Coord {
    i32 qry,ref; //quer, reference
};

//struct Alignment {
//    std::vector<Coord> path;
//    float score;
//
//    static void pybind(py::module_ &m) {
//        c.def_property_readonly("path", [](Alignment &d) -> pybind11::array_t<Coord> {
//             return pybind11::array_t<Trace>(d.path_.size(), d.path_.data());
//        });
//        c.def_readonly("score", &Alignment::score);
//    }
//};

//ReadAln coords_to_aln(const std::vector<Coords> &coords, const ProcessedRead &read) {
//
//}

#ifdef PYBIND
//PyArray<Coord> get_guided_bands(PyArray<i32> bc_refs, PyArray<i32> bc_samps, PyArray<i32> event_samps);
#endif


struct DtwDF : public DataFrame<int, int, float, float> {
    static constexpr NameArray names = {"start", "length", "current", "stdv"}; 
    ColType<0> &start = std::get<0>(data_);  
    ColType<1> &length = std::get<1>(data_);  
    ColType<2> &current = std::get<2>(data_);  
    ColType<3> &stdv = std::get<3>(data_);  
    using DataFrame::DataFrame;                              

    DtwDF(const std::vector<Coord> &path, const ProcessedRead &read, size_t ref_len) : DtwDF(ref_len) {
        i32 i = 0, prev_ref = -1, qry_st = -1;

        //Event evt{-1,-1,-1,-1};
        for (auto itr = path.rbegin(); itr != path.rend(); itr++) {
            if (itr->ref > prev_ref) {
                if (prev_ref >= 0) {
                    auto evt = read.merge_events(qry_st, itr->qry);
                    //std::cout << i << " " << evt.start << " " << evt.mean << " ";
                    start[i] = evt.start;
                    length[i] = evt.length;
                    current[i] = evt.mean;
                    stdv[i] = evt.stdv;
                    i++;
                }
                //std::cout << ((itr->ref)-prev_ref) << "\n"; 
                prev_ref = itr->ref;
                qry_st = itr->qry;
            }
        }
        auto evt = read.merge_events(qry_st, path.front().qry);
        start[i] = evt.start;
        length[i] = evt.length;
        current[i] = evt.mean;
        stdv[i] = evt.stdv;

        //std::cout.flush();
    }

    void set_signal(const ProcessedRead &read) {
        if (current.size() == 0) {
            current = ValArray<float>(height);
            stdv = ValArray<float>(height);
        }

        for (size_t i = 0; i < height; i++) {
            auto seg = static_cast<std::valarray<float>>(read.signal[std::slice(start[i], length[i], 1)]);
            current[i] = seg.sum() / seg.size();
            auto deltas = seg - current[i];
            stdv[i] = sqrt((deltas*deltas).sum() / seg.size());
        }
    }

    static py::class_<DtwDF> pybind(py::module_ &m) {
        auto c = DataFrame::pybind<DtwDF>(m, "_DtwDF");
        c.def("set_signal", &DtwDF::set_signal);
        return c;
    }
};



template <typename ModelType>
class BandedDTW {
    public:

    using KmerType = typename ModelType::kmer_t;

    public:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const DtwParams PRMS;
    const DTWCostFn cost_fn_;

    //const PyArray<float> &qry_vals_;
    const ProcessedRead &qry_vals_;
    size_t event_start_, event_end_;
    const std::vector<KmerType> &ref_vals_;
    const PyArray<Coord> &ll_;

    const ModelType model_;

    //i32 PRMS.band_width;

    //std::vector< std::pair<i32,i32> > ll_;
    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector<Coord> path_;
    float score_sum_;

    static DTWCostFn get_cost_fn(const std::string &cost_str) {
        if (cost_str == "abs_diff") return DTWCostFn::ABS_DIFF;
        if (cost_str == "norm_pdf") return DTWCostFn::NORM_PDF;
        if (cost_str == "z_score") return DTWCostFn::Z_SCORE;
        std::cerr << "Error: unknown DTW cost function: \""
                  << cost_str << "\". Defaulting to \"abs_diff\"\n";
        return DTWCostFn::ABS_DIFF;
    }

    public:
    BandedDTW(const DtwParams &prms,
              //const PyArray<float> &qry_vals,   
              const ProcessedRead &read,
              size_t event_start, size_t event_end,
              const std::vector<KmerType> &ref_vals,
              const ModelType &model,
              const PyArray<Coord> &ll) :
            PRMS(prms),
            cost_fn_(get_cost_fn(PRMS.cost_fn)),
            qry_vals_(read), //(qry_vals),
            event_start_(event_start), event_end_(event_end),
            ref_vals_(ref_vals),
            ll_(ll),
            model_(model) {
        init_mat();
        fill_mat();          
        traceback();
    }

    u32 band_count() const {
        return ref_size() + qry_size() - 1;
    }

    template <typename T>
    bool in_range(i32 i, const PyArray<T> &v) {
        return static_cast<u32>(i) >= 0 && static_cast<u32>(i) < v.size();
    }

    template <typename T>
    bool in_range(i32 i, const std::vector<T> &v) {
        return static_cast<u32>(i) >= 0 && static_cast<u32>(i) < v.size();
    }

    size_t ref_size() const {
        return (ref_vals_.size());
    }

    size_t qry_size() const {
        return (event_end_-event_start_);
    }

    //void init_mat(const std::vector<bool> &rmoves) {
    void init_mat() {
        //Set first two bands to inf
        //Maybe could work around, but might not matter for streaming
        mat_.resize(PRMS.band_width * band_count());
        bcrumbs_.resize(mat_.size());

        //Set origin (0,0) to 0
        auto k0 = 0, k1 = PRMS.band_width;
        for (; k0 < PRMS.band_width; k0++) {
            mat_[k0] = mat_[k1++] = MAX_COST;
        }

        i32 orig = -ll_[0].ref;
        mat_[orig] = 0;
    }

    i32 band_ref_end(u32 band) const {
        return ll_[band].ref + PRMS.band_width;
    }

    i32 band_start(size_t band) const {
        return PRMS.band_width * band;
    }

    u32 band_coord(u32 band, u32 offs) const {
        return band_start(band) + offs;
    }

    struct MatCoords {
        size_t band, k; //Band index, matrix index
        i32 i,j; //Event coord, ref coord
        i32 v, h, adj_min, adj_max,
            d, diag_min, diag_max;
    };

    void fill_mat() {

        for (size_t band = 2; band < band_count(); band++) {
            i32 q = ll_[band].qry,
                r = ll_[band].ref;
            
            u32 band_start = PRMS.band_width * band,
                band_end = band_start + PRMS.band_width;
            
            //adjacent (left/up) band coordinates
            i32 aband  = band - 1, 
                astart = aband * PRMS.band_width,
                aend   = astart + PRMS.band_width,
                ak     = astart + (r - ll_[aband].ref - 1);

            //diagonal band coordinates
            i32 dband  = band - 2, 
                dstart = dband * PRMS.band_width,
                dend   = dstart + PRMS.band_width,
                dk     = dstart + (r - ll_[dband].ref - 1);

            //std::cout << mat_.size() << " "
            //          << q << " " 
            //          << r << " " 
            //          << ak << " " 
            //          << dk << "\n";

            for (size_t k = band_start; k < band_end; k++) {

                //TODO can I compute starting locaiton from ll_?
                if (q >= 0 & q < qry_size() && in_range(r, ref_vals_)) {
                    float cost = costfn(qry_vals_.events[event_start_+q].mean, ref_vals_[r]);
                    float ds,hs,vs;

                    if (ak >= astart && ak <= aend) { 
                        hs = mat_[ak] + cost * PRMS.stay_cost;
                    } else hs = MAX_COST;

                    if (ak+1 >= astart && ak+1 <= aend) {
                        vs = mat_[ak+1] + cost * PRMS.skip_cost;
                    } else vs = MAX_COST;

                    if (dk >= dstart && dk <= dend) {
                        ds = mat_[dk] + cost * PRMS.move_cost;
                    } else ds = MAX_COST;

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


                } else {
                    mat_[k] = MAX_COST;
                }

                q -= 1;
                r += 1;
                ak += 1;
                dk += 1;
            }
        }
    }

    //i32 up_k(i32 k0) {

    //    u32 band = k0
    //}

    //i32 diag_k(u32 k0) {

    //}

    void traceback() {

        i32 path_band = band_count()-1,
            offs = ref_size() - ll_[path_band].ref - 1,
            path_k = band_coord(path_band, offs);
        score_sum_ = mat_[path_k];

        while (path_band > 1) {

            auto path_ll = ll_[path_band];
            auto next_ll = path_ll;

            auto offs = path_k - band_start(path_band);
            //std::cout << path_band << " " << path_k << " " << offs << "\n";
            path_.push_back({path_ll.qry-offs+event_start_, path_ll.ref+offs});

            i32 shift = 0;

            switch(bcrumbs_[path_k]) {
            case Move::D:
                //std::cout << "D";
                path_band -= 2;
                next_ll = ll_[path_band];
                shift = 2 * PRMS.band_width + (next_ll.ref - path_ll.ref + 1);
                break;

            case Move::H:
                //std::cout << "H";
                path_band -= 1;
                next_ll = ll_[path_band];
                shift = PRMS.band_width + (next_ll.ref - path_ll.ref + 1);
                break;
            case Move::V:
                //std::cout << "V";
                path_band -= 1;
                next_ll = ll_[path_band];
                shift = PRMS.band_width + (next_ll.ref - path_ll.ref);
                break;
            }
            //std::cout << path_band << "\t"
            //          << path_k << "\t"
            //          << shift << "\n";
            path_k -= shift;
        }

        auto c = path_.back();
        c.qry = event_start_;
        while (c.ref > 0) {
            c.ref -= 1;
            path_.push_back(c);
        }

        //std::cout << "\n";
        //path_.push_back({event_start_,0});
    }

    std::vector<Coord> get_path() {
        return path_;
    }

    DtwDF get_aln() {
        return DtwDF(path_, qry_vals_, ref_size());
    }

    void fill_aln(Alignment<ModelType> &aln) {
        aln.dtw = AlnDF(aln.seq.mpos);
        auto &dtw = aln.dtw;

        i32 i = 0, prev_ref = -1, qry_st = -1;
        //Event evt{-1,-1,-1,-1};
        for (auto itr = path_.rbegin(); itr != path_.rend(); itr++) {
            if (itr->ref > prev_ref) {
                if (prev_ref >= 0) {
                    auto evt = qry_vals_.merge_events(qry_st, itr->qry);
                    dtw.samples.append(evt.start, evt.start+evt.length);
                    dtw.current[i] = evt.mean;
                    dtw.current_sd[i] = evt.stdv;
                    i++;
                }
                //std::cout << ((itr->ref)-prev_ref) << "\n"; 
                prev_ref = itr->ref;
                qry_st = itr->qry;
            }
        }
        auto evt = qry_vals_.merge_events(qry_st, path_.front().qry);
        dtw.samples.append(evt.start, evt.start+evt.length);
        dtw.current[i] = evt.mean;
        dtw.current_sd[i] = evt.stdv;

    }

    float score() {
        return score_sum_;
    }

    float mean_score() {
        return score_sum_ / path_.size();
    }

    std::vector<float> get_flat_mat() const {
        std::vector<float> flat_mat(qry_size()*ref_size(), 0);

        for (size_t band = 0; band < band_count(); band++) {

            auto q0 = ll_[band].qry,    r0 = ll_[band].ref,
                 q1 = q0 - PRMS.band_width, r1 = r0 + PRMS.band_width; 

            size_t k0 = band_start(band),
                   k1 = k0 + PRMS.band_width;

            u32 head_clip = 0,
                tail_clip = 0;

            if (q0 >= qry_size()) {
                head_clip = q0 - qry_size() + 1;
            }
            if (r0 < 0) {//TODO could be else if?
                head_clip = -r0;
            }

            if (r1 > ref_size()) {
                tail_clip = r1 - ref_size();
            }
            if (q1 < 0) {
                tail_clip = -q1;
            }

            q0 -= head_clip;
            r0 += head_clip;
            k0 += head_clip;
            q1 += tail_clip;
            r1 -= tail_clip;
            k1 -= tail_clip;

            float band_min = MAX_COST, band_max = 0;
            for (size_t k = k0; k < k1; k++) {
                if (mat_[k] < MAX_COST && mat_[k] > band_max) {
                    band_max = mat_[k];
                }
                band_min = std::min(band_min, mat_[k]);
            }
            float band_span = band_max-band_min;

            auto q = q0;
            auto r = r0;
            for (size_t k = k0; k < k1; k++) {
                size_t i = ref_size()*q + r;
                if (mat_[k] < MAX_COST) {
                    flat_mat[i] = 1-((mat_[k]-band_min) / band_span);
                } else {
                    flat_mat[i] = 0;
                }
                q -= 1;
                r += 1;
            }
        }

        return flat_mat;
    }

    protected:

    float costfn(float current, KmerType kmer) {
        switch(cost_fn_) {
            case DTWCostFn::NORM_PDF:
            return -model_.norm_pdf(current, kmer);

            case DTWCostFn::Z_SCORE:
            return model_.z_score(current, kmer);

            default:
            return model_.abs_diff(current, kmer);
        }
    }

    public:

    #ifdef PYBIND
    #define PY_BANDED_DTW_METH(P) c.def(#P, &BandedDTW<ModelType>::P);
    //static void pybind_defs(pybind11::class_<BandedDTW<ModelType>> &c) {
    static void pybind_defs(pybind11::module_ &m, const std::string &suffix) {
        pybind11::class_<BandedDTW<ModelType>> c(m, ("BandedDTW" + suffix).c_str());

        c.def(pybind11::init<const DtwParams &,
                             //const PyArray<float>&, 
                             const ProcessedRead&, 
                             size_t, size_t,
                             const std::vector<KmerType>&,
                             const ModelType&,
                             const PyArray<Coord>&>());

        //c.def(pybind11::init<const DtwParams &, 
        //                     const std::vector<float>&, 
        //                     const std::vector<KmerType>&,
        //                     const ModelType&,
        //                     const std::vector<std::pair<i32,i32>>& >());

        PY_BANDED_DTW_METH(get_path)
        PY_BANDED_DTW_METH(get_aln)
        PY_BANDED_DTW_METH(score)
        PY_BANDED_DTW_METH(mean_score)
        //c.def_readonly("ll", BandedDTW<ModelType>::ll_);
        PY_BANDED_DTW_METH(mean_score)
        PY_BANDED_DTW_METH(get_flat_mat)
        PY_BANDED_DTW_METH(fill_aln)

        c.def_property_readonly("ll", [](BandedDTW<ModelType> &d) -> pybind11::array_t<Coord> {
             return pybind11::array_t<Coord>(d.ll_.size(), d.ll_.data);
        });

        c.def_property_readonly("path", [](BandedDTW<ModelType> &d) -> pybind11::array_t<Coord> {
             return pybind11::array_t<Coord>(d.path_.size(), d.path_.data());
        });
    }
    #endif
};

//template <typename ModelType>
//class StaticBDTW : public BandedDTW<ModelType> {
//
//    public:
//
//    using KmerType = typename ModelType::kmer_t;
//
//    StaticBDTW(
//           const DtwParams &prms,
//           const std::vector<float> &qry_vals,   
//           const std::vector<KmerType> &ref_vals,
//           const ModelType &model) : 
//            BandedDTW<ModelType>(prms, qry_vals, ref_vals, model) {
//
//        init_mat();
//        BandedDTW<ModelType>::fill_mat();
//        BandedDTW<ModelType>::traceback();
//    }
//
//    void init_mat() {
//        auto slope = static_cast<float>(BandedDTW<ModelType>::ref_size()) / BandedDTW<ModelType>::qry_size();
//        auto shift = static_cast<i32>(BandedDTW<ModelType>::PRMS.band_width * BandedDTW<ModelType>::PRMS.band_shift);
//        i32 q = 0, r = 0;
//
//        //std::cout << qry_size() << "\t" << ref_size() << " SIZE\n";
//        //std::cout << slope << "\t" << shift << " SLOPE\n";
//
//        for (size_t i = 0; i < BandedDTW<ModelType>::band_count(); i++) {
//            BandedDTW<ModelType>::ll_.push_back({q+shift,r-shift});
//
//            float down_dist  = std::abs( (r - slope * (q+1)) ),
//                  right_dist = std::abs( ((r+1) - slope * q) );
//            if (down_dist < right_dist) q += 1;
//            else r += 1;
//        }
//
//        //std::cout << ll_.front().qry << "\t"
//        //          << ll_.front().ref << "\n"
//        //          << ll_.back().qry << "\t"
//        //          << ll_.back().ref << "\n";
//
//        BandedDTW<ModelType>::init_mat();
//    }
//
//    #ifdef PYBIND
//    //static void pybind_defs(pybind11::class_<StaticBDTW, BandedDTW<ModelType>> &c) {
//    static void pybind_defs(pybind11::module_ &m, const std::string &suffix) {
//        pybind11::class_<StaticBDTW, BandedDTW<ModelType>> c(m, ("StaticBDTW" + suffix).c_str());
//        c.def(pybind11::init<const DtwParams &,
//                             const std::vector<float>&, 
//                             const std::vector<KmerType>&,
//                             const ModelType &>());
//    }
//    #endif
//};

#endif
