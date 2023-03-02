#ifdef PYBIND

#ifndef ALN_HPP
#define ALN_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "pore_model.hpp"
namespace py = pybind11;

//struct RefCoords : public IntervalIndex<i64> {
//    RefCoords
//}

struct AlnDF {
    IntervalIndex<i64> index;
    IntervalIndex<i32> samples;
    ValArray<float> current, current_sd; 

    AlnDF() {}

    AlnDF(IntervalIndex<i64> index_) : index(index_) {
        //current(index.length),
        //current_sd(index.length)
        current = ValArray<float>(index_.length);
        current_sd = ValArray<float>(index_.length);
    }

    AlnDF(IntervalIndex<i64> index_, py::array_t<i32> start_, py::array_t<i32> length_, py::array_t<float> current_, py::array_t<float> current_sd_) : 
        index(index_),
        samples(start_, length_),
        current(init_arr(current_)),
        current_sd(init_arr(current_sd_)) {}

    AlnDF(IntervalIndex<i64> index_, IntervalIndex<i32> &samples_, py::array_t<float> current_, py::array_t<float> current_sd_) : 
        index(index_),
        samples(samples_),
        current(init_arr(current_)),
        current_sd(init_arr(current_sd_)) {}

    AlnDF slice(size_t i, size_t j) {
        AlnDF ret(index.islice(i, j));
        for (size_t k = i; k < j; k++) {
            ret.samples.append(samples.coords[k]);
            ret.current[k-i] = current[k];
            ret.current_sd[k-i] = current_sd[k];
        }
        return ret;
    }

    bool empty() const {
        return size() == 0;
    }

    size_t size() const {
        return index.length;
    }

    template <typename T>
    static ValArray<T> init_arr(py::array_t<T> &a) {
        auto info = a.request();
        return ValArray<T>(static_cast<T*>(info.ptr), static_cast<size_t>(info.shape[0]));
    }

    static py::class_<AlnDF> pybind(py::module_ &m) {
        py::class_<AlnDF> c(m, "_AlnDF");
        c.def(py::init<IntervalIndex<i64> &>());
        c.def(py::init<IntervalIndex<i64>&, IntervalIndex<i32>&, py::array_t<float>, py::array_t<float>>());
        c.def(py::init<IntervalIndex<i64>&, py::array_t<i32>, py::array_t<i32>, py::array_t<float>, py::array_t<float>>());
        c.def("slice", &AlnDF::slice);
        c.def("__len__", &AlnDF::size);
        c.def_readwrite("index", &AlnDF::index);
        c.def_readwrite("samples", &AlnDF::samples);
        c.def_readwrite("current", &AlnDF::current);
        c.def_readwrite("current_sd", &AlnDF::current_sd);
        return c;
    }
};

struct CmpDF {
    
    IntervalIndex<i64> index;
    ValArray<float> dist, jaccard;

    CmpDF() {}

    CmpDF(AlnDF &a, AlnDF &b) {
        //std::cout << "init\n";
        //std::cout.flush();

        index = a.index.intersect(b.index);

        dist.resize(index.length);
        jaccard.resize(index.length);

        size_t i = 0, j = 0, k = 0;

        int ref;
        //std::cout << "loop\n";
        //std::cout.flush();

        //TODO could optimize so its a single linear scan
        //each get_index call is logn, but n is pretty small
        for (size_t i = 0; i < index.length; i++) {
            //std::cout << "top " << i << "\n";
            //std::cout.flush();
            auto ref = index[i];
            auto ai = a.index.get_index(ref), bi = b.index.get_index(ref);
            auto ac = a.samples.coords[ai];
            auto bc = b.samples.coords[bi];
            //std::cout << "mid " << i << "\n";
            //std::cout.flush();
            if (ac.is_valid() && bc.is_valid()) {
                jaccard[i] = calc_jaccard(ac.start, ac.end, bc.start, bc.end);
            } else {
                jaccard[i] = jaccard.NA;
            }
            //std::cout << "bot " << i << "\n";
            //std::cout.flush();
        }

        //std::cout << "predist\n";
        //std::cout.flush();

        //ref = rec[0].ref;
        DistIter ab(a,b), ba(b,a);

        //for (auto &r : rec) {
        for (i = 0; i < index.length; i++) {
            //std::cout << "top2 " << i << "\n";
            //std::cout.flush();
            auto ref = index[i];
            //std::cout << "next2 " << i << "\n";
            //std::cout.flush();
            auto d0 = ab.next_dist(ref);
            //std::cout << "d0 " << i << "\n";
            //std::cout.flush();
            auto d1 = ba.next_dist(ref);
            //std::cout << "d1 " << i << "\n";
            //std::cout.flush();
            auto denom = d0.length + d1.length;
            //std::cout << r.ref << " (" << d0.dist << ", " << d0.length << ") "
            //          << "(" << d1.dist << ", " << d1.length << ")\n";
            if (denom > 0) {
                dist[i] = (d0.dist + d1.dist) / denom;
            } else {
                dist[i] = dist.NA;
                if (ab.ended() && ba.ended()) break;
            }
            //std::cout << "bot " << i << "\n";
            //std::cout.flush();
        }
    }

    size_t size() const {
        return index.length;
    }

    bool empty() const {
        return index.length == 0;
    }

    struct DistIter {
        AlnDF &a, &b;
        size_t i, j;

        DistIter(AlnDF &a_, AlnDF &b_) : 
            a{a_}, b{b_}, i{0}, j{0} {}

        struct Coefs {float dist, length;};
        Coefs next_dist(int ref) {
            Coefs c = {0,0};
            for (; i < a.index.length; i++) {
                if (a.index[i] == ref) break;
                else if (a.index[i] > ref) {
                    return c;
                }
            }
            auto ac = a.samples.coords[i];
            if (!ac.is_valid()) return c;

            auto bcs = b.samples.coords;

            for (; j < b.index.length; j++) {
                if (bcs[j].is_valid() && bcs[j].end > ac.start) break;
            }
            if (ended() || bcs[j].start >= ac.end) { 
                return c;
            }

            int len,dist;
            do {
                len = std::min(ac.end, bcs[j].end) - std::max(ac.start, bcs[j].start);
                assert(len > 0);
                c.dist += len * std::abs(a.index[i] - b.index[j]);
                c.length += len;
                //std::cout << i << " " << j << " " << c.dist << " " << c.length << "\n";
                j++;
            } while (j < b.index.length && bcs[j].start < ac.end);

            do {
                j--;
            } while (j > 0 && j < b.samples.length && bcs[j].start == bcs[j-1].start);

            return c;
        }

        bool ended() const {
            return i == a.size() || j == b.size();
        }
    };

    float calc_jaccard(int start_a, int end_a, int start_b, int end_b) {
        auto in = std::min(end_a, end_b) - std::max(start_a, start_b);
        if (in < 0) return 1;
        auto un = static_cast<float>(std::max(end_a, end_b) - std::min(start_a, start_b));
        return 1 - (in / un);
    }

    static void pybind(py::module_ &m) {
        auto c = py::class_<CmpDF>(m, "_CmpDF");
        c.def(py::init<AlnDF&, AlnDF&>());
        c.def("empty", &CmpDF::empty);
        c.def("__len__", &CmpDF::size);
        c.def_readonly("index", &CmpDF::index);
        c.def_readonly("dist", &CmpDF::dist);
        c.def_readonly("jaccard", &CmpDF::jaccard);
    }
};

template <typename ModelType>
struct Alignment {
    std::string read_id;
    Sequence<ModelType> seq;
    AlnDF dtw, moves;
    CmpDF mvcmp;

    Alignment(const std::string &read_id_, Sequence<ModelType> seq_) :
        read_id(read_id_), seq(seq_) {
    }

    void set_dtw(AlnDF df) {
        dtw = df;
    }

    void set_moves(AlnDF df) {
        moves = df;
    }

    void calc_mvcmp() {
        if (dtw.empty() || moves.empty()) {
            throw std::runtime_error("Both 'dtw' and 'moves' must be defined");
        }
        mvcmp = CmpDF(dtw, moves);
    }

    static void pybind(py::module_ &m, std::string suffix) {
        py::class_<Alignment> c(m, ("_Alignment"+suffix).c_str());
        c.def(py::init<const std::string &, Sequence<ModelType>>());
        c.def("set_dtw", &Alignment::set_dtw);
        c.def("set_moves", &Alignment::set_moves);
        c.def("calc_mvcmp", &Alignment::calc_mvcmp);
        c.def_readonly("read_id", &Alignment::read_id);
        c.def_readonly("seq", &Alignment::seq);
        c.def_readonly("_mvcmp", &Alignment::mvcmp);
        c.def_readonly("_dtw", &Alignment::dtw);
        c.def_readonly("_moves", &Alignment::moves);
    }
};


#endif
#endif
