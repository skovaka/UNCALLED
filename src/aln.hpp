#ifdef PYBIND

#ifndef ALN_HPP
#define ALN_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "seq.hpp"
#include "event_detector.hpp"
namespace py = pybind11;

struct AlnDF {
    IntervalIndex<i64> index;
    IntervalIndex<i32> samples;
    IntervalIndex<i32> sample_bounds;
    ValArray<float> current, current_sd, events; 

    AlnDF() {}

    AlnDF(IntervalIndex<i64> index_) : index(index_) {
        current = ValArray<float>(index_.length);
        current_sd = ValArray<float>(index_.length);
    }

    AlnDF(IntervalIndex<i64> index_, py::array_t<i32> start_, py::array_t<i32> length_, py::array_t<float> current_, py::array_t<float> current_sd_) : 
        index(index_),
        samples(start_, length_),
        current(init_arr(current_)),
        current_sd(init_arr(current_sd_)) {
        init_sample_bounds();
    }

    AlnDF(IntervalIndex<i64> index_, IntervalIndex<i32> &samples_, py::array_t<float> current_, py::array_t<float> current_sd_) : 
        index(index_),
        samples(samples_),
        current(init_arr(current_)),
        current_sd(init_arr(current_sd_)) {
        init_sample_bounds();
    }
    
    void init_sample_bounds() {
        sample_bounds = IntervalIndex<i32>();
        sample_bounds.append(samples.get_start(), samples.get_end());
    }

    void set_hard_gaps(const std::vector<size_t> idxs) {
        sample_bounds = {};
        auto start = samples.get_start(), end = start;
        for (auto i : idxs) {
            end = samples.coords[i-1].end;
            assert(start < end);
            sample_bounds.append(start, end);
            start = samples.coords[i].start;
        }
        end = samples.get_end();
        sample_bounds.append(start, end);
    }

    AlnDF slice(size_t i, size_t j) {
        //std::cout << "in\n";
        //std::cout.flush();
        AlnDF ret(index.islice(i, j));
        //std::cout << "it\n";
        //std::cout.flush();
        for (size_t k = i; k < j; k++) {
            ret.samples.append(samples.coords[k]);
            ret.current[k-i] = current[k];
            ret.current_sd[k-i] = current_sd[k];
        }
        //std::cout << sample_bounds.to_string() << "\n";
        //std::cout << "SL" << ret.samples.get_start() << " " << ret.samples.get_end() << "\n";
        //std::cout.flush();
        ret.sample_bounds = sample_bounds.slice(ret.samples.get_start(), ret.samples.get_end());
        //std::cout << "rt\n";
        //std::cout.flush();
        return ret;
    }

    AlnDF trim_na() {
        size_t start_na = 0, end_na = 0;
        for (size_t i = 0; i < size(); i++) {
            if (index.coords[i].is_valid()) break;
            start_na += 1;
        }
        for (size_t i = size()-1; i < size(); i--) {
            if (index.coords[i].is_valid()) break;
            end_na += 1;
        }
        return slice(start_na, size()-end_na);
    }

    void set_signal(const ProcessedRead &read) {
        if (current.size() == 0) {
            current = ValArray<float>(size());
            current_sd = ValArray<float>(size());
        }

        for (size_t i = 0; i < size(); i++) {
            auto c = samples.coords[i];
            auto seg = static_cast<std::valarray<float>>(read.signal[std::slice(c.start, c.length(), 1)]);
            current[i] = seg.sum() / seg.size();
            auto deltas = seg - current[i];
            current_sd[i] = sqrt((deltas*deltas).sum() / seg.size());
        }
    }


    void mask_i(size_t i) {
        if (current.size() > 0) current[i] = current.NA;         
        if (current_sd.size() > 0) current_sd[i] = current_sd.NA;
        samples.coords[i].clear();                               
    }

    void mask(py::array_t<bool> mask_py) {
        PyArray<bool> mask(mask_py);
        if (mask.size() != size()) {
            throw std::runtime_error("Mask must be same length as dataframe");
        }
        for (size_t i = 0; i < size(); i++) {
            if (!mask[i]) {
                mask_i(i);
                //if (current.size() > 0) current[i] = current.NA;
                //if (current_sd.size() > 0) current_sd[i] = current_sd.NA;
                //samples.coords[i].clear();
            }
        }
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
        c.def("trim_na", &AlnDF::trim_na);
        c.def("empty", &AlnDF::empty);
        c.def("mask", &AlnDF::mask);
        c.def("set_signal", &AlnDF::set_signal);
        c.def("__len__", &AlnDF::size);
        c.def_readwrite("index", &AlnDF::index);
        c.def_readwrite("samples", &AlnDF::samples);
        c.def_readwrite("current", &AlnDF::current);
        c.def_readwrite("current_sd", &AlnDF::current_sd);
        c.def_readwrite("events", &AlnDF::events);
        c.def_readwrite("sample_bounds", &AlnDF::sample_bounds);
        return c;
    }
};

struct CmpDF {
    
    IntervalIndex<i64> index;
    ValArray<float> dist, jaccard, current, current_sd;

    CmpDF() {}

    CmpDF(AlnDF &a, AlnDF &b, bool symetric) {

        index = a.index.intersect(b.index);

        dist.resize(index.length);
        jaccard.resize(index.length);

        if (!a.current.empty() && !b.current.empty()) {
            current.resize(index.length);
        }
        if (!a.current_sd.empty() && !b.current_sd.empty()) {
            current_sd.resize(index.length);
        }

        i64 ref;

        //TODO could optimize so its a single linear scan
        //each get_index call is logn, but n is pretty small
        for (size_t i = 0; i < index.length; i++) {
            auto ref = index[i];
            auto ai = a.index.get_index(ref), bi = b.index.get_index(ref);
            auto ac = a.samples.coords[ai];
            auto bc = b.samples.coords[bi];
            if (ac.is_valid() && bc.is_valid()) {
                jaccard[i] = calc_jaccard(ac.start, ac.end, bc.start, bc.end);
            } else {
                jaccard[i] = jaccard.NA;
            }
            if (!a.current.empty() && !b.current.empty()) {
                current[i] = a.current[ai] - b.current[bi];
            }
            if (!a.current_sd.empty() && !b.current_sd.empty()) {
                current_sd[i] = a.current_sd[ai] - b.current_sd[bi];
            }
        }

        if (symetric) {
            DistIter ab(a,b), ba(b,a);

            //for (auto &r : rec) {
            for (size_t i = 0; i < index.length; i++) {
                auto ref = index[i];
                auto d0 = ab.next_dist(ref);
                auto d1 = ba.next_dist(ref);
                auto denom = d0.length + d1.length;
                if (denom > 0) {
                    dist[i] = (d0.dist + d1.dist) / denom;
                } else {
                    dist[i] = dist.NA;
                    if (ab.ended() && ba.ended()) break;
                }
            }
        } else {
            DistIter ab(a,b);

            //for (auto &r : rec) {
            for (size_t i = 0; i < index.length; i++) {
                auto ref = index[i];
                auto d = ab.next_dist(ref);
                if (d.length > 0) {
                    dist[i] = d.dist / d.length;
                } else {
                    dist[i] = dist.NA;
                    if (ab.ended()) break;
                }
            }

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
        c.def(py::init<AlnDF&, AlnDF&, bool>());
        c.def("empty", &CmpDF::empty);
        c.def("__len__", &CmpDF::size);
        c.def_readonly("index", &CmpDF::index);
        c.def_readonly("dist", &CmpDF::dist);
        c.def_readonly("jaccard", &CmpDF::jaccard);
        c.def_readonly("current", &CmpDF::current);
        c.def_readonly("current_sd", &CmpDF::current_sd);
    }
};

template <typename ModelType>
struct Alignment {
    int id;
    std::string read_id;
    Sequence<ModelType> seq;
    AlnDF dtw, moves;
    CmpDF mvcmp, dtwcmp;

    Alignment(int id_, const std::string &read_id_, Sequence<ModelType> seq_) :
        id(id_), read_id(read_id_), seq(seq_) {
    }

    void set_dtw(AlnDF df) {
        dtw = df;
    }

    void set_moves(AlnDF df) {
        moves = df;
    }

    void calc_dtwcmp(Alignment &aln) {
        if (!(dtw.empty() || aln.dtw.empty())) {
            dtwcmp = CmpDF(dtw, aln.dtw, true);
        }
    }

    void calc_mvcmp() {
        if (!(dtw.empty() || moves.empty())) {
            mvcmp = CmpDF(dtw, moves, false);
            //mvcmp = CmpDF(dtw, moves, true);
        }
    }

    void mask_skips(bool keep_best) {
        if (keep_best) mask_skips<true>();
        else mask_skips<false>();
    }

    void mask_skips(py::array_t<float> cost_py) {
        PyArray<float> costs(cost_py);
        
        auto prev_start = dtw.samples.coords[0].start;
        int nskips = 0;
        for (size_t i = 1; i < dtw.size(); i++) {
            auto start = dtw.samples.coords[i].start;
            if (start == prev_start) {
                nskips++;
            } else if (nskips > 0) {
                mask_skip(i-nskips-1, i, costs);
                nskips = 0;
            }
            prev_start = start;
        }

        if (nskips > 0) {
            mask_skip(dtw.size()-nskips-1, dtw.size(), costs);
        }
    }

    void mask_skip(size_t start, size_t end, const PyArray<float> &cost) {
        float best = INFINITY;
        size_t best_i = -1;
        for (size_t i = start; i < end; i++) {
            if (cost[i] < best) {
                best_i = i;
                best = cost[i];
            }
        }
        for (size_t i = start; i < end; i++) {
            if (i != best_i) {
                dtw.mask_i(i);
            }
        }
    }

    template<bool keep_best>
    void mask_skips() {
        auto prev_start = dtw.samples.coords[0].start;
        int nskips = 0;
        for (size_t i = 1; i < dtw.size(); i++) {
            auto start = dtw.samples.coords[i].start;
            if (start == prev_start) {
                nskips++;
            } else if (nskips > 0) {
                mask_skip<keep_best>(i-nskips-1, i);
                nskips = 0;
            }
            prev_start = start;
        }

        if (nskips > 0) {
            mask_skip<keep_best>(dtw.size()-nskips-1, dtw.size());
        }
    }

    template <bool keep_best> 
    typename std::enable_if<keep_best, void>::type
    mask_skip(size_t start, size_t end) {
        float best = INFINITY;
        size_t best_i = -1;
        for (size_t i = start; i < end; i++) {
            auto abs_diff = abs(dtw.current[i] - seq.current[i]);
            if (abs_diff < best) {
                best_i = i;
                best = abs_diff;
            }
        }
        for (size_t i = start; i < end; i++) {
            if (i != best_i) {// && abs(dtw.current[i] - seq.current[i]) > 2*dtw.current_sd[i]) {
            dtw.mask_i(i);
            }
        }
    }

    template <bool keep_best> 
    typename std::enable_if<!keep_best, void>::type
    mask_skip(size_t start, size_t end) {
        for (size_t i = start; i < end; i++) {
            auto abs_diff = abs(dtw.current[i] - seq.current[i]);
            //if (abs_diff >= seq.current_sd[i]) {
            if (abs_diff >= 2*dtw.current_sd[i]) {
                dtw.mask_i(i);
            }
        }
    }


    static void pybind(py::module_ &m, std::string suffix) {
        py::class_<Alignment> c(m, ("_Alignment"+suffix).c_str());
        c.def(py::init<int, const std::string &, Sequence<ModelType>>());
        c.def("set_dtw", &Alignment::set_dtw);
        c.def("set_moves", &Alignment::set_moves);
        c.def("calc_dtwcmp", &Alignment::calc_dtwcmp);
        c.def("calc_mvcmp", &Alignment::calc_mvcmp);
        c.def("mask_skips", static_cast<void (Alignment::*)(bool)>(&Alignment::mask_skips));
        c.def("mask_skips", static_cast<void (Alignment::*)(py::array_t<float>)>(&Alignment::mask_skips));
        c.def_readwrite("id", &Alignment::id);
        c.def_readwrite("read_id", &Alignment::read_id);
        c.def_readonly("seq", &Alignment::seq);
        c.def_readonly("_dtwcmp", &Alignment::dtwcmp);
        c.def_readonly("_mvcmp", &Alignment::mvcmp);
        c.def_readonly("_dtw", &Alignment::dtw);
        c.def_readonly("_moves", &Alignment::moves);
    }
};


#endif
#endif
