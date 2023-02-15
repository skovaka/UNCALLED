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

PYBIND11_MAKE_OPAQUE(std::valarray<float>);
PYBIND11_MAKE_OPAQUE(std::valarray<i8>);
PYBIND11_MAKE_OPAQUE(std::valarray<i16>);
PYBIND11_MAKE_OPAQUE(std::valarray<i32>);
PYBIND11_MAKE_OPAQUE(std::valarray<i64>);
PYBIND11_MAKE_OPAQUE(std::valarray<u8>);
PYBIND11_MAKE_OPAQUE(std::valarray<u16>);
PYBIND11_MAKE_OPAQUE(std::valarray<u32>);
PYBIND11_MAKE_OPAQUE(std::valarray<u64>);

template<typename T>
void pybind_valarray(py::module_ &m, std::string suffix) {
    auto name = "Valarray"+suffix;
    py::class_<std::valarray<T>> a(m, name.c_str(), py::buffer_protocol());
    a.def_buffer([](std::valarray<T> &c) -> py::buffer_info { 
        return py::buffer_info( 
            &c[0], 
            sizeof(T), 
            py::format_descriptor<T>::format(), 
            1, 
            {c.size()}, 
            {sizeof(T)} 
        ); 
    });
    a.def("__repr__", [](std::valarray<T> &a) -> std::string { 
        std::stringstream ss;
        ss << "[";
        if (a.size() > 6) {
            for (size_t i = 0; i < 3; i++) {
                ss << a[i] << " ";
            }
            ss << "...";
            for (size_t i = a.size()-4; i < a.size(); i++) {
                ss << " " << a[i];
            }
        } else {
            for (size_t i = 0; i < a.size()-1; i++) {
                ss << a[i] << " ";
            }
            ss << a[a.size()-1];
        }
        ss << "]";
        return ss.str();
    });

    a.def("__getitem__", static_cast< T& (std::valarray<T>::*)(size_t)> (&std::valarray<T>::operator[]));
    a.def("__len__", &std::valarray<T>::size);
    a.def("__iter__", [](const std::valarray<T> &v) { 
            return py::make_iterator(std::begin(v), std::end(v)); 
        }, py::keep_alive<0, 1>())
    ;


    a.def("to_numpy", [](std::valarray<T> &a) -> py::array_t<T> { 
        return py::array_t<T>{static_cast<py::ssize_t>(a.size()), &a[0]}; 
    });
}

template<typename T>
struct PyArray {
    py::buffer_info info;
    T *data;
    size_t size_;

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

template <typename... Types>
class DataFrame {
    public:
    static constexpr size_t width = sizeof...(Types);

    using DataTuple = std::tuple<std::valarray<Types>...>;
    using NameArray = std::array<const char *, width>;
    DataTuple data_;

    template <size_t I>
    using ColType = typename std::tuple_element<I, DataTuple>::type;

    const size_t height;

    DataFrame(const DataFrame &) = delete;
    DataFrame &operator=(const DataFrame &) = delete;
    DataFrame(DataFrame&&) = default;

    DataFrame(py::array_t<Types>... arrays) : 
        data_  { init_arr(arrays)... },
        height {init_height()} {
    }

    DataFrame(size_t length) : 
        data_  { std::valarray<Types>(length)... },
        height {init_height()} {
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
    static std::valarray<T> init_arr(py::array_t<T> a) {
        auto info = a.request();
        return std::valarray<T>(static_cast<T*>(info.ptr), static_cast<size_t>(info.shape[0]));
    }

    public:


    template<class Subclass>
    static py::class_<Subclass> pybind(py::module_ &m, const char *name) {
        py::class_<Subclass> c(m, name);
        c.def(py::init<py::array_t<Types>...>());
        c.def(py::init<size_t>());

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

template <typename T>
struct Interval {
    T start, end;

    Interval(T s, T e) : start(s), end(e) {}
    Interval(const std::pair<T,T>& p) : Interval(p.first, p.second) {}

    static const T NA = std::numeric_limits<T>::max();

    bool contains(T val) const {
        return val >= start && val < end;
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
class IntervalIndex {
    public:

    std::vector<Interval<T>> coords;
    std::vector<size_t> starts;
    size_t length;

    static constexpr size_t MAX_IDX = -1;//std::numeric_limits<size_t>::max();

    IntervalIndex(std::vector<std::pair<T,T>> coords_, bool sorted=false) 
            : coords(coords_.begin(), coords_.end()) {

        if (!sorted) {
            std::sort(coords.begin(), coords.end());
        }

        starts.reserve(coords.size());
        length = 0;
        for (auto c : coords) {
            starts.push_back(length);
            length += c.length();
        }
    }

    T operator[] (size_t i) {
        if (i > length) throw std::out_of_range("Interval index of range");
        size_t j = std::upper_bound(starts.begin(), starts.end(), i) - starts.begin() - 1;
        return coords[j].start + (i - starts[j]);
    }

    std::string to_string() const {
        std::stringstream ss;
        for (auto c : coords) {
            ss << c.to_string() << " ";
        }
        return ss.str();
    }

    Interval<T> get_interval(size_t i) const {
        return coords[i];
    }

    std::valarray<T> get_lengths() const {
        std::valarray<T> ret(size());
        for (size_t i = 0; i < size(); i++) {
            ret[i] = coords[i].end - coords[i].start;
        }
        return ret;
    }

    std::valarray<T> expand() const {
        std::valarray<T> ret(length);
        size_t i = 0;
        for (auto &intv : coords) {
            for (T v = intv.start; v < intv.end; v++) {
                ret[i++] = v;
            }
        }
        return ret;
    }

    size_t get_index(T val) const {
        Interval<T> q = {val, Interval<T>::NA};
        auto itr = std::lower_bound(coords.begin(), coords.end(), q);

        if (itr > coords.begin()) {
            size_t i = static_cast<size_t>(itr - coords.begin()) - 1;
            if (coords[i].contains(val)) {
                return starts[i] + (val - coords[i].start);
            }
        }
        return MAX_IDX;
    }

    size_t size() const {
        return coords.size();
    }

    #ifdef PYBIND
    static void pybind(py::module_ &m, const std::string &suffix) {
        PYBIND11_NUMPY_DTYPE(Interval<T>, start, end);

        py::class_<IntervalIndex>(m, ("IntervalIndex"+suffix).c_str())
            .def(py::init<std::vector<std::pair<T,T>>>())
            .def("__getitem__", py::vectorize(&IntervalIndex::operator[]))
            .def("__len__", &IntervalIndex::size)
            .def("__repr__", &IntervalIndex<T>::to_string)
            .def("expand", &IntervalIndex::expand)
            .def("get_interval", py::vectorize(&IntervalIndex::get_interval))
            .def("get_index", py::vectorize(&IntervalIndex::get_index))
            .def_property_readonly("lengths", &IntervalIndex::get_lengths)
            ;
            
        py::class_<Interval<T>>(m, ("Interval"+suffix).c_str())
            .def("__repr__", &Interval<T>::to_string)
            .def("__getitem__", &Interval<T>::operator[])
            .def_readwrite("start", &Interval<T>::start)
            .def_readwrite("end", &Interval<T>::end);
    }
    #endif
};

#endif
#endif
