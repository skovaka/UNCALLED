#ifdef PYBIND

#ifndef DATAFRAME_HPP
#define DATAFRAME_HPP

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

template<typename T>
struct PyArray {

    std::vector<T> data_vec;
    py::buffer_info info;
    T *data;
    size_t size_;

    PyArray(size_t length, T fill) :
        data_vec { length, fill },
        info { data_vec.data(), length },
        data { data_vec.data() },
        size_ { length } {}

    PyArray(T *ptr, size_t length) :
        info { ptr, length },
        data { ptr },
        size_ { length } {}

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

    using DataTuple = std::tuple<PyArray<Types>...>;
    using NameArray = std::array<const char *, width>;
    DataTuple data_;

    template <size_t I>
    using ColType = typename std::tuple_element<I, DataTuple>::type;

    const size_t height;

    DataFrame(const DataFrame &) = delete;
    DataFrame &operator=(const DataFrame &) = delete;
    DataFrame(DataFrame&&) = default;

    DataFrame(py::array_t<Types>... arrays) : 
        data_  { PyArray<Types>(arrays)... },
        height {init_height()} {
    }

    size_t init_height() {
        auto size = std::get<0>(data_).size();
        return init_height<1>(size);
    }

    template <size_t I>
    typename std::enable_if<I < width, size_t>::type
    init_height(size_t size) {
        if (std::get<I>(data_).size() != size) {
            throw std::runtime_error("All DataFrame columns must be same size");
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


    template<class Subclass>
    static void pybind(py::module_ &m, const char *name) {
        py::class_<Subclass> c(m, name);
        c.def(py::init<py::array_t<Types>...>());

        c.def("__len__", [](Subclass &c) -> size_t {return c.height;});
        c.attr("names") = py::cast(Subclass::names);
        c.attr("width") = py::cast(Subclass::width);
        c.def_readonly("height", &Subclass::height);

        pybind_col<Subclass>(c);
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

#endif
#endif
