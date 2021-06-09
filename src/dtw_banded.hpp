#ifndef _INCL_DTW_BANDED
#define _INCL_DTW_BANDED

#include <vector>
#include <iostream>
#include <cfloat>
#include "util.hpp"
#include "pore_model.hpp"

//TODO don't duplicate klen
//probably should just set constant, or declare elsewhere
const KmerLen DTWB_KLEN = KmerLen::k5;

class BandedDTW {
    public:

    struct Coord {
        i32 qry,ref; //quer, reference
    };

    //using 

    //protected:
    public:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const std::vector<float> qry_vals_;
    const std::vector<u16> ref_vals_;

    const PoreModel<DTWB_KLEN> model_;

    i32 band_width_;

    //std::vector< std::pair<i32,i32> > ll_;
    std::vector<float> mat_;
    std::vector<Coord> ll_;
    std::vector<Move> bcrumbs_;
    std::vector<Coord> path_;
    float score_sum_;

    public:
    BandedDTW(const std::vector<float> &qry_vals,   
              const std::vector<u16> &ref_vals,
              const PoreModel<DTWB_KLEN> &model,
              u32 band_width) :
            qry_vals_(qry_vals),
            ref_vals_(ref_vals),
            model_(model),
            band_width_(band_width) {}

    BandedDTW(const std::vector<float> &qry_vals,   
              const std::vector<u16> &ref_vals,
              const PoreModel<DTWB_KLEN> &model,
              u32 band_width,
              const std::vector<Coord> &ll) :
                BandedDTW(qry_vals, ref_vals, model, band_width) {

        ll_ = ll;

        init_mat();
        fill_mat();          
        traceback();
    }

    BandedDTW(const std::vector<float> &qry_vals,   
              const std::vector<u16> &ref_vals,
              const PoreModel<DTWB_KLEN> &model,
              u32 band_width,
              const std::vector<std::pair<i32,i32>> &ll) :
                BandedDTW(qry_vals, ref_vals, model, band_width) {
        //assert(lower.size() == left.size());

        ll_.resize(ll.size());
        //for (size_t i = 0; i < lower.size(); i++) {
        for (size_t i = 0; i < ll.size(); i++) {
            ll_[i] = {ll[i].first, ll[i].second};
        }

        init_mat();
        fill_mat();          
        traceback();
    }


    u32 band_count() const {
        return ref_size() + qry_size() - 1;
    }

    template <typename T>
    bool in_range(i32 i, const std::vector<T> &v) {
        return static_cast<u32>(i) >= 0 && static_cast<u32>(i) < v.size();
    }

    i32 ref_size() const {
        return static_cast<i32>(ref_vals_.size());
    }

    i32 qry_size() const {
        return static_cast<i32>(qry_vals_.size());
    }

    //void init_mat(const std::vector<bool> &rmoves) {
    void init_mat() {
        //Set first two bands to inf
        //Maybe could work around, but might not matter for streaming
        mat_.resize(band_width_ * band_count());
        bcrumbs_.resize(mat_.size());

        //Set origin (0,0) to 0
        auto k0 = 0, k1 = band_width_;
        for (; k0 < band_width_; k0++) {
            mat_[k0] = mat_[k1++] = MAX_COST;
        }

        i32 orig = -ll_[0].ref;
        mat_[orig] = 0;
    }

    i32 band_ref_end(u32 band) const {
        return ll_[band].ref + band_width_;
    }

    i32 band_start(size_t band) const {
        return band_width_ * band;
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
            
            u32 band_start = band_width_ * band,
                band_end = band_start + band_width_;
            
            //adjacent (left/up) band coordinates
            i32 aband  = band - 1, 
                astart = aband * band_width_,
                aend   = astart + band_width_,
                ak     = astart + (r - ll_[aband].ref - 1);

            //diagonal band coordinates
            i32 dband  = band - 2, 
                dstart = dband * band_width_,
                dend   = dstart + band_width_,
                dk     = dstart + (r - ll_[dband].ref - 1);

            for (size_t k = band_start; k < band_end; k++) {

                //TODO can I compute starting locaiton from ll_?
                if (in_range(q, qry_vals_) && in_range(r, ref_vals_)) {
                    float cost = costfn(qry_vals_[q], ref_vals_[r]);
                    float ds,hs,vs;

                    if (ak >= astart && ak <= aend) { 
                        hs = mat_[ak] + cost;
                    } else hs = MAX_COST;

                    if (ak+1 >= astart && ak+1 <= aend) {
                        vs = mat_[ak+1] + cost;
                    } else vs = MAX_COST;

                    if (dk >= dstart && dk <= dend) {
                        ds = mat_[dk] + cost;
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
            path_.push_back({path_ll.qry-offs, path_ll.ref+offs});

            i32 shift = 0;

            switch(bcrumbs_[path_k]) {
            case Move::D:
                //std::cout << "D";
                path_band -= 2;
                next_ll = ll_[path_band];
                shift = 2 * band_width_ + (next_ll.ref - path_ll.ref + 1);
                break;

            case Move::H:
                //std::cout << "H";
                path_band -= 1;
                next_ll = ll_[path_band];
                shift = band_width_ + (next_ll.ref - path_ll.ref + 1);
                break;
            case Move::V:
                //std::cout << "V";
                path_band -= 1;
                next_ll = ll_[path_band];
                shift = band_width_ + (next_ll.ref - path_ll.ref);
                break;
            }
            //std::cout << path_band << "\t"
            //          << path_k << "\t"
            //          << shift << "\n";
            path_k -= shift;
        }

        //std::cout << "\n";
        path_.push_back({0,0});
    }

    std::vector<Coord> get_path() {
        return path_;
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
                 q1 = q0 - band_width_, r1 = r0 + band_width_; 

            size_t k0 = band_start(band),
                   k1 = k0 + band_width_;

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

    float costfn(float pA, u16 kmer) {
        return abs(pA-model_.get_mean(kmer));
        //return -model_.match_prob(pA,kmer);
    }

    public:

    #ifdef PYBIND
    #define PY_BANDED_DTW_METH(P) c.def(#P, &BandedDTW::P);
    static void pybind_defs(pybind11::class_<BandedDTW> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTWB_KLEN>&,
                             u32,
                             const std::vector<Coord>&>());

        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTWB_KLEN>&,
                             u32,
                             const std::vector<std::pair<i32,i32>>& >());
        PY_BANDED_DTW_METH(get_path)
        PY_BANDED_DTW_METH(score)
        PY_BANDED_DTW_METH(mean_score)
        c.def_readonly("ll", &BandedDTW::ll_);
        PY_BANDED_DTW_METH(mean_score)
        PY_BANDED_DTW_METH(get_flat_mat)

        c.def_property_readonly("ll", [](BandedDTW &d) -> pybind11::array_t<Coord> {
             return pybind11::array_t<Coord>(d.ll_.size(), d.ll_.data());
        });

        c.def_property_readonly("path", [](BandedDTW &d) -> pybind11::array_t<Coord> {
             return pybind11::array_t<Coord>(d.path_.size(), d.path_.data());
        });
        PYBIND11_NUMPY_DTYPE(Coord, qry, ref);
    }
    #endif
};

class StaticBDTW : public BandedDTW {

    public:

    float band_center_;

    StaticBDTW(const std::vector<float> &qry_vals,   
           const std::vector<u16> &ref_vals,
           const PoreModel<DTWB_KLEN> &model,
           u32 band_width,
           float band_center=0.5) : 
            BandedDTW(qry_vals, ref_vals, model, band_width) {
        band_center_ = band_center;

        init_mat();
        fill_mat();
        traceback();
    }

    void init_mat() {
        auto slope = static_cast<float>(ref_size()) / qry_size();
        auto shift = static_cast<i32>(band_width_ * band_center_);
        i32 q = 0, r = 0;

        //std::cout << qry_size() << "\t" << ref_size() << " SIZE\n";
        //std::cout << slope << "\t" << shift << " SLOPE\n";

        for (size_t i = 0; i < band_count(); i++) {
            ll_.push_back({q+shift,r-shift});

            float down_dist  = std::abs( (r - slope * (q+1)) ),
                  right_dist = std::abs( ((r+1) - slope * q) );
            if (down_dist < right_dist) q += 1;
            else r += 1;
        }

        //std::cout << ll_.front().qry << "\t"
        //          << ll_.front().ref << "\n"
        //          << ll_.back().qry << "\t"
        //          << ll_.back().ref << "\n";

        BandedDTW::init_mat();
    }

    #ifdef PYBIND
    static void pybind_defs(pybind11::class_<StaticBDTW, BandedDTW> &c) {
        c.def(pybind11::init<const std::vector<float>&, 
                             const std::vector<u16>&,
                             const PoreModel<DTWB_KLEN> &,
                             u32, float>());
    }
    #endif
};

#endif
