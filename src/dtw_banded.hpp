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
    struct Trace {
        u64 qry, ref;
    };

    //protected:
    public:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const PoreModel<DTWB_KLEN> model_;

    const std::vector<u16> ref_vals_;
    const std::vector<float> qry_vals_;

    i32 band_width_, band_count_;

    std::vector< std::pair<i32,i32> > ll_;
    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector<Trace> path_;
    float score_sum_;

    public:
    BandedDTW(const std::vector<float> &col_vals,   
         const std::vector<u16> &row_vals,
         const std::vector<bool> &rmoves,
         u32 band_width,
         const PoreModel<DTWB_KLEN> &model) :
             model_(model),
             ref_vals_(row_vals),
             qry_vals_(col_vals),
             band_width_(band_width),
             //band_count_(rmoves.size()) {
             band_count_(ref_len() + qry_len()) {

        mat_.resize(band_width_ * band_count_);
        bcrumbs_.resize(mat_.size());

        //std::cout << qry_len() << "\t" << ref_len() << " LEN\n";

        Timer t;

        init_mat(rmoves);
        //std::cout << "Init: " << t.lap() << "\n";
        compute_matrix();          
        //std::cout << "Fill: " << t.lap() << "\n";
        traceback();
        //std::cout << "Trace: " << t.lap() << "\n";

    }

    template <typename T>
    bool in_range(i32 i, std::vector<T> v) {
        return static_cast<u32>(i) >= 0 && static_cast<u32>(i) < v.size();
    }

    i32 ref_len() const {
        return static_cast<i32>(ref_vals_.size());
    }

    i32 qry_len() const {
        return static_cast<i32>(qry_vals_.size());
    }

    void init_mat(const std::vector<bool> &rmoves) {
        u32 shift = band_width_/2;
        i32 i = shift, j = -shift;
        for (auto move_right : rmoves) {
            //std::cout << ll_.size() << "\t" 
            //          << (i-band_width_) << "," << (j+band_width_) << "\t" << i << "," << j << "\n";
            ll_.push_back({i,j});
            i += not move_right;
            j += move_right;
        }

        //Set first two bands to inf
        //Maybe could work around, but might not matter for streaming
        auto k0 = 0, k1 = band_width_;
        for (; k0 < band_width_; k0++) {
            mat_[k0] = mat_[k1++] = MAX_COST;
        }

        //Set origin (0,0) to 0
        i32 orig = -ll_[0].second;
        mat_[orig] = 0;
    }

    i32 band_ref_end(u32 band) {
        return ll_[band].second + band_width_;
    }

    u32 band_start(size_t band) {
        return band_width_ * band;
    }

    u32 band_coord(u32 band, u32 offs) {
        return band_start(band) + offs;
    }

    struct MatCoords {
        size_t band, k; //Band index, matrix index
        i32 i,j; //Event coord, ref coord
        i32 v, h, adj_min, adj_max,
            d, diag_min, diag_max;
    };

    void compute_matrix() {

        for (size_t band = 2; band < band_count_; band++) {
            i32 i = ll_[band].first,
                j = ll_[band].second;

            u32 band_start = band * band_width_,
                band_end = band_start + band_width_;
            
            //adjacent (left/up) band coordinates
            i32 aband  = band - 1, 
                astart = aband * band_width_,
                aend   = astart + band_width_,
                ak     = astart + (j - ll_[aband].second - 1);

            //diagonal band coordinates
            i32 dband  = band - 2, 
                dstart = dband * band_width_,
                dend   = dstart + band_width_,
                dk     = dstart + (j - ll_[dband].second - 1);

            for (size_t k = band_start; k < band_end; k++) {

                //TODO can I compute starting locaiton from ll_?
                if (in_range(i, qry_vals_) && in_range(j, ref_vals_)) {
                    float cost = costfn(qry_vals_[i], ref_vals_[j]);
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

                i -= 1;
                j += 1;
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
        score_sum_ = MAX_COST;
        i32 path_band = band_count_, path_k = band_width_;
        for (size_t band = band_count_-1; band > 0; band--) {
    
            //Find index with last ref coordinate
            //stop search if not in band
            auto offs = ref_len() - ll_[band].second - 1;
            //std::cout << offs << " AH\n";
            if (offs >= band_width_) break;

            auto k = band_coord(band, offs);
            //std::cout << ll_[band].first-offs << " " << ll_[band].second+offs << " " << mat_[k] << " AH\n";

            if (mat_[k] < score_sum_) {
                path_band = band;
                path_k = k;
                score_sum_ = mat_[k];
            }
        }


        while (path_band > 1) {

            auto path_ll = ll_[path_band];
            auto next_ll = path_ll;

            auto offs = path_k - band_start(path_band);
            //std::cout << path_band << " " << path_k << " " << offs << "\n";
            path_.push_back({path_ll.first-offs, path_ll.second+offs});

            i32 shift = 0;

            switch(bcrumbs_[path_k]) {
            case Move::D:
                path_band -= 2;
                next_ll = ll_[path_band];
                shift = 2 * band_width_ + (next_ll.second - path_ll.second + 1);
                break;

            case Move::H:
                path_band -= 1;
                next_ll = ll_[path_band];
                shift = band_width_ + (next_ll.second - path_ll.second + 1);
                break;
            case Move::V:
                path_band -= 1;
                next_ll = ll_[path_band];
                shift = band_width_ + (next_ll.second - path_ll.second);
                break;
            }
            //std::cout << path_band << "\t"
            //          << path_k << "\t"
            //          << shift << "\n";
            path_k -= shift;
        }

        path_.push_back({0,0});
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
                             const std::vector<bool>&,
                             u32,
                             const PoreModel<DTWB_KLEN> &>());
        PY_BANDED_DTW_METH(get_path)
        PY_BANDED_DTW_METH(score)
        PY_BANDED_DTW_METH(mean_score)
        c.def_readonly("ll", &BandedDTW::ll_);
        PY_BANDED_DTW_METH(mean_score)

        c.def_property_readonly("path", [](BandedDTW &d) -> pybind11::array_t<Trace> {
             return pybind11::array_t<Trace>(d.path_.size(), d.path_.data());
        });
        PYBIND11_NUMPY_DTYPE(Trace, qry, ref);
    }
    #endif
};

#endif
