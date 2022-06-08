#ifdef PYBIND

#ifndef COMPARE_HPP
#define COMPARE_HPP

#include <limits>
#include <cmath>

#include "dataframe.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

const auto INF = std::nanf("");

using AlnCoords = PyArray<AlnCoord>;

class Compare {
    public:
    struct Rec {
        int ref;
        float jaccard, mean_ref_dist;
    };

    using CoordPair = std::pair<size_t, size_t>;

    std::vector<Rec> rec;

    Compare(PyArray<AlnCoord> &a, PyArray<AlnCoord> &b) {
        size_t i = 0, j = 0;

        int ref;

        while (i < a.size && j < b.size) {
            float jaccard = INF;

            if (i == a.size) {
                ref = b[j++].ref;
            } else if (j == b.size) {
                ref = a[i++].ref;
            } else if (a[i].ref < b[j].ref) {
                ref = a[i++].ref;     
            } else if (a[i].ref > b[j].ref) {
                ref = b[j++].ref;
            } else {
                ref = a[i].ref;
                jaccard = calc_jaccard(a[i].start, a[i].end, b[j].start, b[j].end);
                i++;j++;

                rec.push_back({ref, jaccard, INF});
            }
        }

        ref = rec[0].ref;
        DistIter ab(a,b), ba(b,a);

        for (auto &r : rec) {
            auto d0 = ab.next_dist(r.ref),
                 d1 = ba.next_dist(r.ref);
            auto denom = d0.length + d1.length;
            //std::cout << r.ref << " (" << d0.dist << ", " << d0.length << ") "
            //          << "(" << d1.dist << ", " << d1.length << ")\n";
            if (denom > 0) {
                r.mean_ref_dist = (d0.dist + d1.dist) / denom;
            } else if (ab.ended() && ba.ended()) break;
        }
    }

    struct DistIter {
        AlnCoords &a, &b;
        size_t i, j;

        DistIter(AlnCoords &a_, AlnCoords &b_) : 
            a{a_}, b{b_}, i{0}, j{0} {}

        struct Coefs {float dist, length;};
        Coefs next_dist(int ref) {
            Coefs c = {0,0};
            for (; i < a.size; i++) {
                if (a[i].ref == ref) break;
                else if (a[i].ref > ref) {
                    //std::cout << "Skipped\n";
                    return c;
                }
            }
            for (; j < b.size; j++) {
                if (b[j].end > a[i].start) break;
            }
            if (ended() || b[j].start >= a[i].end) { 
                //std::cout << "Ended? " << ended() << "\n";
                return c;
            }

            int len,dist;
            do {
                len = std::min(a[i].end, b[j].end) - std::max(a[i].start, b[j].start);
                assert(len > 0);
                c.dist += len * std::abs(a[i].ref - b[j].ref);
                c.length += len;
                //std::cout << i << " " << j << " " << c.dist << " " << c.length << "\n";
                j++;
            } while (j < b.size && b[j].start < a[i].end);

            do {
                j--;
            } while (j > 0 && j < b.size && b[j].start == b[j-1].start);

            return c;
        }

        bool ended() const {
            return i == a.size || j == b.size;
        }
    };

    //def weighted_dist(idxs, idx_b, swap):
    //   coord_a = merge.loc[idxs]
    //    coord_b = merge.loc[idx_b]
    //    samp_overlaps = coord_a["end_"+a].clip(upper=coord_b["end_"+b]) - coord_a["start_"+a].clip(lower=coord_b["start_"+b])
    //    ref_dists = np.abs(idxs.get_level_values(0).to_numpy()-idx_b[0])
    //    return np.sum(samp_overlaps), np.sum(ref_dists * samp_overlaps)

    float calc_jaccard(int start_a, int end_a, int start_b, int end_b) {
        auto in = std::min(end_a, end_b) - std::max(start_a, start_b);
        if (in < 0) return 1;
        auto un = static_cast<float>(std::max(end_a, end_b) - std::min(start_a, start_b));
        return 1 - (in / un);
    }

    static void pybind(py::module_ &m) {
        PYBIND11_NUMPY_DTYPE(Rec, ref, jaccard, mean_ref_dist);
        auto c = py::class_<Compare>(m, "Compare");
        c.def(py::init<AlnCoords&, AlnCoords&>());
        c.def("to_numpy", [](Compare &c) -> py::array_t<Rec> {
            return py::array_t<Rec>{c.rec.size(), c.rec.data()};
        });
    }

};

#endif
#endif
