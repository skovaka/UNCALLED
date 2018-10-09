#ifndef DTW_HPP
#define DTW_HPP

#include <vector>
#include <iostream>
#include <cfloat>
#include "basepairs.hpp"

enum SubSeqDTW {NONE, ROW, COL};

template < class RowT, class ColT, typename Func >
class DTW {
    private:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const std::vector<RowT> &rvals_;
    const std::vector<ColT> &cvals_;

    const Func &C;

    SubSeqDTW subseq_;
    float weights_[3];

    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector< std::pair<u64, u64> > path_;
    float score_sum_;

    public:
    DTW(const std::vector<RowT> &row_vals, 
        const std::vector<ColT> &col_vals,
        const Func &cost_fn,
        //std::function<float(RowT, ColT)> cost_fn,
        SubSeqDTW which_subseq = SubSeqDTW::NONE,
        float d_wt=1, float h_wt=1, float v_wt=1) :

          rvals_(row_vals),
          cvals_(col_vals),
          C(cost_fn),
          subseq_(which_subseq),
          weights_{d_wt, h_wt, v_wt},
          mat_(std::vector<float>(row_vals.size() * col_vals.size())),
          bcrumbs_(std::vector<Move>(mat_.size())) {
    
        compute_matrix();          
        traceback();
    }


    inline float dscore(u64 i, u64 j) {
        if (j > 0 && i > 0) return mat_[cvals_.size()*(i-1) + j-1];
        else if ((j == i) || 
                 (i == 0 && subseq_ == SubSeqDTW::COL) || 
                 (j == 0 && subseq_ == SubSeqDTW::ROW)) return 0;
        return MAX_COST;
    }

    inline float hscore(u64 i, u64 j) {
        if (j > 0) return mat_[cvals_.size()*i + j-1];
        else if (subseq_ == SubSeqDTW::ROW) return 0;
        return MAX_COST;
    }

    inline float vscore(u64 i, u64 j) {
        if (i > 0) return mat_[cvals_.size()*(i-1) + j];
        else if (subseq_ == SubSeqDTW::COL) return 0;
        else return MAX_COST;
    }

    void compute_matrix() {
        u64 k = 0;
        float cost, ds, hs, vs;
        for (u64 i = 0; i < rvals_.size(); i++) {
            for (u64 j = 0; j < cvals_.size(); j++) {

                cost = C(rvals_[i], cvals_[j]);
                ds = dscore(i,j) + cost * weights_[Move::D];
                hs = hscore(i,j) + cost * weights_[Move::H];
                vs = vscore(i,j) + cost * weights_[Move::V];
                //std::cout << cost << " " << ds << " " << hs << " " << vs << std::endl;

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

                //std::cout << i << " " << j << " " << mat_[k-1] << std::endl;
            }
        }
    }

    void traceback() {

        u64 i = rvals_.size()-1, j = cvals_.size()-1;

        switch (subseq_) {
            case SubSeqDTW::ROW:
                for (u64 k = 0; k < rvals_.size(); k++) {
                    //std::cout << k << " " << mat_[k*cvals_.size() + j] << "\n";
                    if (mat_[k*cvals_.size() + j] < mat_[i*cvals_.size() + j]) { 
                        i = k;
                    }
                }
                break;
            case SubSeqDTW::COL:
                for (u64 k = 0; k < cvals_.size(); k++) {
                    //std::cout << k << " " << mat_[i*cvals_.size() + k] << "\n";
                    if (mat_[i*cvals_.size() + k] < mat_[i*cvals_.size() + j]) {
                        j = k;
                    }
                }
                break;
            default:
                break;
        }

        score_sum_ = mat_[i*cvals_.size() + j];

        path_.push_back(std::pair<u64, u64>(i, j));

        u64 k = i*cvals_.size() + j;
        while(!(i == 0 || subseq_ == SubSeqDTW::ROW) ||
              !(j == 0 || subseq_ == SubSeqDTW::COL)) {

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

            path_.push_back(std::pair<u64, u64>(i, j));
        }

    }

    float score() {
        return score_sum_;
    }

    float mean_score() {
        return score_sum_ / path_.size();
    }

    //friend std::ostream &operator<< (std::ostream &out, const DTW &a);
    void print_path() {
        for (auto p = path_.rbegin(); p != path_.rend(); p++) {
            std::cout << p->first << "\t"
                      << p->second << "\t"
                      << rvals_[p->first] << "\t"
                      << cvals_[p->second]
                      << std::endl;
        }
    } 
};




#endif
