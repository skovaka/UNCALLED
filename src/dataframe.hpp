#ifdef PYBIND

#ifndef DATAFRAME_HPP
#define DATAFRAME_HPP

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "util.hpp"
namespace py = pybind11;

template<typename T>
class ValArray: public std::valarray<T> {
    public:
    using std::valarray<T>::valarray;
    using super = std::valarray<T>;

    static constexpr T init_na() {
        return  std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    }

    static constexpr T NA = init_na();//std::numeric_limits<T>::max();

    T &at(size_t i) {
        return super::operator[](i);
    }

    std::string to_string() {
        std::stringstream ss;
        ss << "[";
        if (super::size() > 6) {
            for (size_t i = 0; i < 3; i++) {
                ss << at(i) << " ";
            }
            ss << "...";
            for (size_t i = super::size()-4; i < super::size(); i++) {
                ss << " " << at(i);
            }
        } else {
            for (size_t i = 0; i < super::size()-1; i++) {
                ss << at(i) << " ";
            }
            ss << at(super::size()-1);
        }
        ss << "]";
        return ss.str();
    }

    py::array_t<T> to_numpy() {
        return py::array_t<T>(static_cast<py::ssize_t>(super::size()), &((*this)[0])); 
    }

    static void pybind(py::module_ &m, std::string suffix) {
        auto name = "ValArray"+suffix;
        py::class_<ValArray> a(m, name.c_str(), py::buffer_protocol());
        a.def_buffer([](ValArray &c) -> py::buffer_info { 
            return py::buffer_info( 
                &c[0], 
                sizeof(T), 
                py::format_descriptor<T>::format(), 
                1, 
                {c.size()}, 
                {sizeof(T)} 
            ); 
        });
        //a.def("__getitem__", static_cast<T& (super::*)(size_t)> (&super::operator[]));
        a.def("__getitem__", &ValArray::at);
        a.def("__repr__", &ValArray::to_string);
        a.def("to_numpy", &ValArray::to_numpy);
        a.def("__len__", &super::size);
        a.def("__iter__", [](const ValArray &v) { 
                return py::make_iterator(std::begin(v), std::end(v)); 
            }, py::keep_alive<0, 1>());

    }
};

template<typename T>
struct PyArray {

    py::buffer_info info;
    T *data;
    size_t size_;

    //TODO turn into ArrayRef, could reference slices
    //PyArray(T *ptr, size_t length) :
    //    info { ptr, length },
    //    data { ptr },
    //    size_ { length } {
    //}

    PyArray(py::array_t<T> arr) :
        info { arr.request() },
        data { static_cast<T*>(info.ptr) },
        size_ { static_cast<size_t>(info.shape[0]) } {}

    PyArray(const PyArray &) = delete;
    PyArray &operator=(const PyArray &) = delete;
    PyArray(PyArray&&) = default;

    T &operator[](size_t i) const {
        return data[i];
    }

    T *begin() const {
        return data;
    }

    T *end() const {
        return &data[size_];
    }

    py::array_t<T> to_numpy() {
        return py::array_t<T>(size_, data);
    }

    size_t size() const {
        return size_;
    }

    template < typename Subclass = PyArray<T> >
    static py::class_<Subclass> pybind(py::module_ &m, const char *name) {
        py::class_<Subclass> c(m, name, py::buffer_protocol());
		c.def(py::init<py::array_t<T>>());
		c.def_buffer([](Subclass &c) -> py::buffer_info {
			return py::buffer_info(
				c.data,                            
				sizeof(T),                          
				py::format_descriptor<T>::format(), 
				1,                                  
				{c.size_},
				{sizeof(T)}
			);
		});
		c.def("to_numpy", &Subclass::to_numpy);

		c.def("__getitem__", &Subclass::operator[]);
		c.def("__len__", &Subclass::size);

        return c;
    }
};


template <typename T>
struct Interval {
    T start, end;
    static const T NA = std::numeric_limits<T>::max();

    Interval() : start(NA), end(NA) {}
    Interval(T s, T e) : start(s), end(e) {}
    Interval(const std::pair<T,T>& p) : Interval(p.first, p.second) {
    }

    void clear() {
        start = NA;
        end = NA;
    }

    bool is_valid() const {
        return !(start == NA || end == NA) && start < end;
    }

    std::pair<T,T> to_pair() const {
        return {start, end};
    }

    bool contains(T val) const {
        return val >= start && val < end;
    }

    Interval intersect(const Interval &b) const {
        return {std::max(start, b.start), std::min(end, b.end)};
    }

    size_t length() const {
        return static_cast<size_t>(end-start);
    }

    T &operator[](size_t i) {
        switch (i) {
            case 0: return start;
            case 1: return end;
        }
        throw std::out_of_range("Interval index out of bounds");
    }

    std::string to_string() const {
        std::stringstream ss;
        ss << "[" << start << ", " << end << ")";
        return ss.str();
    }
};

template <typename T>
bool operator <(const Interval<T> &a, const Interval<T> &b) { 
    return a.start < b.start || (a.start == b.start && a.end < b.end);
}

template <typename T>
bool operator ==(const Interval<T> &a, const Interval<T> &b) { 
    return a.start == b.start && a.end == b.end;
}

template <typename T>
bool operator !=(const Interval<T> &a, const Interval<T> &b) { 
    return a.start != b.start || a.end != b.end;
}

template <typename T>
class IntervalIndex {
    public:

    std::vector<Interval<T>> coords;
    std::vector<size_t> starts;
    size_t length = 0;

    static constexpr size_t MAX_IDX = -1;//std::numeric_limits<size_t>::max();

    IntervalIndex() {}

    IntervalIndex(std::vector<std::pair<T,T>> coords_, bool sorted=false)
            : coords(coords_.begin(), coords_.end()) {

        if (!sorted) {
            std::sort(coords.begin(), coords.end());
        }

        starts.reserve(coords.size());
        length = 0;
        for (auto c : coords) {
            if (c.is_valid()) {
                starts.push_back(length);
                length += c.length();
            }
        }
    }

    IntervalIndex(py::array_t<i32> starts_py, py::array_t<i32> lengths_py) {
        PyArray<i32> starts_(starts_py), lengths_(lengths_py);
        if (starts_.size() != lengths_.size()) {
            throw std::runtime_error("Interval arrays must be same length");
        }

        coords.reserve(starts_.size());
        starts.reserve(starts_.size());
        length = 0;
        for (size_t i = 0; i < starts_.size(); i++) {
            append(starts_[i], starts_[i]+lengths_[i]);
        }
    }

    void init_coords(std::vector<std::pair<T,T>> coords_) {
        //TODO
    }

    void shift(i64 delta) {
        for (auto &c : coords) {
            if (c.is_valid()) {
                c.start += delta;
                c.end += delta;
            }
        }
    }

    IntervalIndex<T> intersect(IntervalIndex<T> other) const {
        IntervalIndex<T> ret;
        auto a = coords.begin(),
             b = other.coords.begin();
        
        while (a != coords.end() && b != other.coords.end()) {
             //std::cout << (a-coords.begin()) << " "
             //         << (b - other.coords.begin()) << "\n";
             if      (b == other.coords.end() || a->end <= b->start) a++;
             else if (a == coords.end() || b->end <= a->start) b++;
             else {
                auto c = a->intersect(*b);
                if (c.is_valid()) {
                    //std::cout << "gotit  " << (a->to_string()) << " " 
                    //          << (b->to_string()) << " "
                    //          << c.to_string() << "\n";
                    ret.append(c);
                } else {
                    std::cout << "failed " << (a->to_string()) << " " 
                              << (b->to_string()) << " "
                              << c.to_string() << "\n";
                }

                if (a->start < b->start) a++;
                else b++;
             }
        }
        return ret;
    }

    IntervalIndex<T> slice(T i, T j) const {
        auto st = get_index(i), en = get_index(j-1);
        return islice(st, en+1);
    }

    IntervalIndex<T> islice(size_t i, size_t j) const {
        IntervalIndex<T> ret;
        auto st = idx_to_interval(i); 
        auto en = idx_to_interval(j-1);

        ret.reserve(en-st+1);

        auto intv = coords[st];
        intv.start += i - starts[st];


        if (st == en) {
            intv.end = intv.start + (j-i);
        } else {
            ret.append(intv);
            for (size_t i = st+1; i < en; i++) {
                ret.append(coords[i]);
            }
            intv = coords[en];
            intv.end = intv.start + (j-i) - ret.length;
        }
        ret.append(intv);

        return ret;
    }

    void reserve(size_t size) {
        coords.reserve(size);
        starts.reserve(size);
    }

    void append(Interval<T> intv) {
        coords.push_back(intv);
        starts.push_back(length);
        if (intv.is_valid()) {
            length += intv.length();
        }
    }

    void append(T start, T end) {
        append({start,end});
    }

    size_t idx_to_interval(size_t i) const {
        if (i > length) throw std::out_of_range("Interval index of range");
        return std::upper_bound(starts.begin(), starts.end(), i) - starts.begin() - 1;
    }

    T operator[] (size_t i) {
        auto j = idx_to_interval(i);
        return coords[j].start + (i - starts[j]);
    }

    std::string to_string() const {
        std::stringstream ss;
        for (auto c : coords) {
            ss << c.to_string() << " ";
        }
        return ss.str();
    }

    ValArray<T> get_gaps() const {
        ValArray<T> ret(coords.size()-1);
        for (size_t i = 0; i < ret.size(); i++) {
            ret[i] = coords[i+1].start - coords[i].end;
        }
        return ret;
    }

    ValArray<T> get_lengths() const {
        ValArray<T> ret(coords.size());
        for (size_t i = 0; i < coords.size(); i++) {
            ret[i] = coords[i].end - coords[i].start;
        }
        return ret;
    }

    ValArray<T> get_lengths_dedup() const { 
        ValArray<T> ret(coords.size());
        auto prev = coords[0];
        ret[0] = prev.length();
        for (size_t i = 1; i < coords.size(); i++) {
            if (coords[i] != prev) {
                ret[i] = coords[i].length();
                prev = coords[i];
            } else {
                ret[i] = 0;
            }
        }
        return ret;
    }

    ValArray<T> get_starts() const {
        ValArray<T> ret(coords.size());
        for (size_t i = 0; i < coords.size(); i++) {
            ret[i] = coords[i].start;
        }
        return ret;
    }

    T get_start() const {
        return coords.front().start;
    }

    T get_end() const {
        return coords.back().end;
    }

    ValArray<T> get_ends() const {
        ValArray<T> ret(coords.size());
        for (size_t i = 0; i < coords.size(); i++) {
            ret[i] = coords[i].end;
        }
        return ret;
    }

    ValArray<T> expand() const {
        ValArray<T> ret(length);
        size_t i = 0;
        for (auto &intv : coords) {
            for (T v = intv.start; v < intv.end; v++) {
                ret[i++] = v;
            }
        }
        return ret;
    }

    size_t get_interval_idx(T val) const {
        Interval<T> q = {val, Interval<T>::NA};
        auto itr = std::lower_bound(coords.begin(), coords.end(), q);

        if (itr > coords.begin()) {
            size_t i = static_cast<size_t>(itr - coords.begin()) - 1;
            if (coords[i].contains(val)) {
                return i;
            }
        }
        return -1;
    }

    Interval<T> get_interval(size_t i) const {
        return coords[i];
    }

    size_t get_index(T val) const {
        auto i = get_interval_idx(val);
        if (i < coords.size()) {
            return starts[i] + (val - coords[i].start);
        } else {
            return MAX_IDX;
        }
    }

    size_t interval_count() const {
        return coords.size();
    }

    #ifdef PYBIND
    static void pybind(py::module_ &m, const std::string &suffix) {
        PYBIND11_NUMPY_DTYPE(Interval<T>, start, end);

        py::class_<IntervalIndex>(m, ("IntervalIndex"+suffix).c_str())
            .def(py::init<std::vector<std::pair<T,T>>>())
            .def("__getitem__", py::vectorize(&IntervalIndex::operator[]))
            //.def("__len__", &IntervalIndex::size)
            .def("interval_count", &IntervalIndex::interval_count)
            .def_readonly("length", &IntervalIndex::length)//TODO bind to __len__
            .def("__repr__", &IntervalIndex<T>::to_string)
            .def("expand", &IntervalIndex::expand)
            .def("shift", &IntervalIndex::shift)
            .def("slice", &IntervalIndex::slice)
            .def("islice", &IntervalIndex::islice)
            .def("intersect", &IntervalIndex::intersect)
            .def("get_interval", py::vectorize(&IntervalIndex::get_interval))
            .def("get_index", py::vectorize(&IntervalIndex::get_index))
            .def_property_readonly("lengths", &IntervalIndex::get_lengths)
            .def_property_readonly("lengths_dedup", &IntervalIndex::get_lengths_dedup)
            .def_property_readonly("gaps", &IntervalIndex::get_gaps)
            .def_property_readonly("starts", &IntervalIndex::get_starts)
            .def_property_readonly("ends", &IntervalIndex::get_ends)
            .def_property_readonly("start", &IntervalIndex::get_start)
            .def_property_readonly("end", &IntervalIndex::get_end)
            ;
            
        py::class_<Interval<T>>(m, ("Interval"+suffix).c_str())
            .def("__repr__", &Interval<T>::to_string)
            .def("__getitem__", &Interval<T>::operator[])
            .def("intersect", &Interval<T>::intersect)
            .def_readwrite("start", &Interval<T>::start)
            .def_readwrite("end", &Interval<T>::end)
            .def("to_tuple", &Interval<T>::to_pair);
    }
    #endif
};

template <typename... Types>
class DataFrame {
    public:
    static constexpr size_t width = sizeof...(Types);

    using DataTuple = std::tuple<ValArray<Types>...>;
    using NameArray = std::array<const char *, width>;
    DataTuple data_;

    template <size_t I>
    using ColType = typename std::tuple_element<I, DataTuple>::type;

    const size_t height;

    DataFrame(const DataFrame &) = default;
    DataFrame &operator=(const DataFrame &) = default;
    DataFrame(DataFrame&&) = default;

    DataFrame(py::array_t<Types>... arrays) : 
        data_  { init_arr(arrays)... },
        height {init_height()} {
    }

    DataFrame(size_t length) : 
        data_  { ValArray<Types>(length)... },
        height { length } {
    }

    DataFrame(size_t length, bool full...) : 
        data_  { init_arr<Types>(length, full)... },
        height { length } {
    }

    size_t init_height() {
        auto size = std::get<0>(data_).size();
        return init_height<1>(size);
    }

    template <size_t I>
    typename std::enable_if<I < width, size_t>::type
    init_height(size_t size) {
        auto arr_size = std::get<I>(data_).size();
        if (arr_size > 0) {
            if (size == 0) {
                size = arr_size;
            } else if (arr_size != size) {
                throw std::runtime_error("All DataFrame columns must be same size or empty");
            }
        }
        return init_height<I+1>(size);
    }

    template <size_t I>
    typename std::enable_if<I == width, size_t>::type
    init_height(size_t size) {
        return size;
    }

    //using 

    template<size_t I, typename T = ColType<I>>
    T &get() {
        return std::get<I>(data_);
    }

    private:

    template <typename T>
    static ValArray<T> init_arr(py::array_t<T> a) {
        auto info = a.request();
        return ValArray<T>(static_cast<T*>(info.ptr), static_cast<size_t>(info.shape[0]));
    }

    template <typename T>
    static ValArray<T> init_arr(size_t n, bool full=true) {
        return full ? ValArray<T>(n) : ValArray<T>();
    }

    public:


    template<class Subclass>
    static py::class_<Subclass> pybind(py::module_ &m, const char *name, bool constructors=true) {
        py::class_<Subclass> c(m, name);

        if (constructors) {
            c.def(py::init<py::array_t<Types>...>());
            c.def(py::init<size_t>());
        }

        c.def("__len__", [](Subclass &c) -> size_t {return c.height;});
        c.attr("names") = py::cast(Subclass::names);
        c.attr("width") = py::cast(Subclass::width);
        c.def_readonly("height", &Subclass::height);

        pybind_col<Subclass>(c);
        
        return c;
    }

    template<class Subclass, size_t I=0, typename T = typename std::tuple_element<I, DataTuple>::type>
    static typename std::enable_if<I < width, void>::type
    pybind_col(py::class_<Subclass> &c) {
        c.def_property_readonly(
            Subclass::names[I], 
            [](Subclass &c) -> T& {
                return std::get<I>(c.data_);
        });
        pybind_col<Subclass, I + 1>(c);
    }

    template<class Subclass, size_t I>
    static typename std::enable_if<I == width, void>::type
    pybind_col(py::class_<Subclass> &c) {
        return;
    }
};

struct AlnCoordsDF : public DataFrame<int, int, int> {
    static constexpr NameArray names = {"ref", "start", "end"}; 
    ColType<0> &ref = std::get<0>(data_);                      
    ColType<1> &start = std::get<1>(data_);                      
    ColType<2> &end = std::get<2>(data_);                      
    using DataFrame::DataFrame;                              
};                

void pybind_dataframes(py::module_ &m);

//TODO 
// redo dataframes so it's
//create ReadAln, contains DtwDF, BcalnDF, and aln info:
//store coordinates (list of ref bound pairs), use to index DFs
//aln_id, read_id, ref_name, ref_start, ref_end, ref_fwd, samp_start, samp_en
//Eventually construct AlnTrack in C++ by concating ReadAlns
//   contains alignments and layers DFs, probably distinct from ReadAln
//   sortable and indexable by aln_id, ref

template<typename T>
struct RecArray : public PyArray<T> {
    using PyArray<T>::PyArray;
    template <typename Subclass>
    static void pybind_rec(py::module_ &m) {
        auto c = RecArray<T>::template pybind<Subclass>(m, Subclass::name);
        c.attr("columns") = py::cast(Subclass::columns);
        Subclass::bind_rec();
    }
};

struct AlnCoord {
    int ref, start, end;
};
struct AlnCoords : public RecArray<AlnCoord> {
    using RecArray<AlnCoord>::RecArray;

    static constexpr const char *name = "_AlnCoords";
    static constexpr std::array<const char *, 3> columns = {"ref", "start", "end"};

    static void bind_rec() {
        PYBIND11_NUMPY_DTYPE(AlnCoord, ref, start, end);
    }
};


#endif
#endif
