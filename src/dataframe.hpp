#ifdef PYBIND
#ifndef DATAFRAME_HPP

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

template<typename T>
struct PyArray {

    py::buffer_info info;
    T *data;
    size_t size;

    PyArray(py::array_t<T> arr) :
        info { arr.request() },
        data { static_cast<T*>(info.ptr) },
        size { static_cast<size_t>(info.shape[0]) } {}

    PyArray(const PyArray &) = delete;
    PyArray &operator=(const PyArray &) = delete;
    PyArray(PyArray&&) = default;

    T &operator[](size_t i) {
        return data[i];
    }

    py::array_t<T> to_numpy() {
        return py::array_t<T>(size, data);
    }

    static void pybind(py::module_ &m, const std::string &suffix) {
        py::class_<PyArray> c(m, ("PyArray"+suffix).c_str(), py::buffer_protocol());
		c.def(py::init<py::array_t<T>>());
		c.def_buffer([](PyArray &c) -> py::buffer_info {
			return py::buffer_info(
				c.data,                            
				sizeof(T),                          
				py::format_descriptor<T>::format(), 
				1,                                  
				{c.size},
				{sizeof(T)}
			);
		});
		c.def("to_numpy", &PyArray::to_numpy);

		c.def("__getitem__", &PyArray::operator[]);
		c.def("__len__", [](PyArray &c) -> size_t {return c.size;});
    }

    private:
};


//TODO maybe just used templatized PyArray + inheritance to process them, cast each dataframe seperately

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
        auto size = std::get<0>(data_).size;
        return init_height<1>(size);
    }

    template <size_t I>
    typename std::enable_if<I < width, size_t>::type
    init_height(size_t size) {
        if (std::get<I>(data_).size != size) {
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

struct AlnCoords : public DataFrame<int, int, int> {
    static constexpr NameArray names = {"ref", "start", "length"}; 
    ColType<0> &a = std::get<0>(data_);                      
    ColType<1> &b = std::get<1>(data_);                      
    ColType<2> &c = std::get<2>(data_);                      
    using DataFrame::DataFrame;                              
};                                                           

void pybind_dataframes(py::module_ &m);


#endif
#endif
