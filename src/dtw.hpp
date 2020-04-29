#ifndef DTW_HPP
#define DTW_HPP

#include <vector>
#include <iostream>
#include <cfloat>
#include "util.hpp"

enum DTWSubSeq {NONE, ROW, COL};


template < class RowT, class ColT, typename Func >
class DTW {
    public:

    typedef struct {
        DTWSubSeq subseq;
        float dw, hw, vw;
        const Func &fn;
    } Prms;


    DTW(const std::vector<RowT> &row_vals,
        const std::vector<ColT> &col_vals,   
        const Prms &p) : 
        prms(p),
        rvals_(row_vals),
        cvals_(col_vals) {

        mat_.resize(rvals_.size() * cvals_.size());
        bcrumbs_.resize(mat_.size());

        compute_matrix();          
        traceback();
    }

    DTW(const std::vector<RowT> &row_vals,
        const std::vector<ColT> &col_vals,   
        DTWSubSeq subseq,
        float dw, 
        float hw, 
        float vw,
        const Func &fn) :
        DTW ( row_vals, col_vals, {subseq,dw,hw,vw,fn} ) {}


    void compute_matrix() {
        u64 k = 0;
        float cost, ds, hs, vs;
        for (u64 i = 0; i < rvals_.size(); i++) {
            for (u64 j = 0; j < cvals_.size(); j++) {

                cost = prms.fn(rvals_[i], cvals_[j]);
                ds = dscore(i,j) + (prms.dw * cost);
                hs = hscore(i,j) + (prms.hw * cost);
                vs = vscore(i,j) + (prms.vw * cost);

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

        switch (prms.subseq) {
            case DTWSubSeq::ROW:
                for (u64 k = 0; k < rvals_.size(); k++) {
                    //std::cout << k << " " << mat_[k*cvals_.size() + j] << "\n";
                    if (mat_[k*cvals_.size() + j] < mat_[i*cvals_.size() + j]) { 
                        i = k;
                    }
                }
                break;
            case DTWSubSeq::COL:
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
        while(!(i == 0 || prms.subseq == DTWSubSeq::ROW) ||
              !(j == 0 || prms.subseq == DTWSubSeq::COL)) {

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
    void print_path() {
        for (auto p = path_.rbegin(); p != path_.rend(); p++) {
            std::cout << p->first << "\t"
                      << p->second << "\t"
                      << rvals_[p->first] << "\t"
                      << cvals_[p->second] << "\t"
                      << "\n";
        }
    } 

    private:
    enum Move {D, H, V}; //Horizontal, vertical, diagonal
    static constexpr float MAX_COST = FLT_MAX / 2.0;

    const Prms prms;

    const std::vector<RowT> &rvals_;
    const std::vector<ColT> &cvals_;

    std::vector<float> mat_;
    std::vector<Move> bcrumbs_;
    std::vector< std::pair<u64, u64> > path_;
    float score_sum_;

    inline float hscore(u64 i, u64 j) {
        if (j > 0) return mat_[cvals_.size()*i + j-1];
        else if (prms.subseq == DTWSubSeq::ROW) return 0;
        return MAX_COST;
    }

    inline float vscore(u64 i, u64 j) {
        if (i > 0) return mat_[cvals_.size()*(i-1) + j];
        else if (prms.subseq == DTWSubSeq::COL) return 0;
        else return MAX_COST;
    }

    inline float dscore(u64 i, u64 j) {
        if (j > 0 && i > 0) return mat_[cvals_.size()*(i-1) + j-1];
        else if ((j == i) || 
                 (i == 0 && prms.subseq == DTWSubSeq::COL) || 
                 (j == 0 && prms.subseq == DTWSubSeq::ROW)) return 0;
        return MAX_COST;
    }
};


#endif
