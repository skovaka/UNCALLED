/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _INCL_UTIL
#define _INCL_UTIL

#include <string>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <cassert>
#include <vector>

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

//Based on github.com/dnbaker/bonsai/blob/master/bonsai/include/util.h
using i8  = std::int8_t;  using u8  = std::uint8_t;
using i16 = std::int16_t; using u16 = std::uint16_t;
using i32 = std::int32_t; using u32 = std::uint32_t;
using i64 = std::int64_t; using u64 = std::uint64_t;

class Timer {
    private:
        std::chrono::high_resolution_clock::time_point start;

    public:
        inline Timer() {
            reset();
        }

        inline void reset() {	
            start = std::chrono::high_resolution_clock::now();
        }

        inline double get() {
            return (std::chrono::duration_cast< std::chrono::duration<double> > (std::chrono::high_resolution_clock::now() - start).count()) * 1000.0;
        }

        inline double lap() {
            double ret = get();
            reset();
            return ret;
        }
};

#define BASE_COUNT 4

const char BASE_CHARS[] {'A', 'C', 'G', 'T', 'N'};

const u8 BASE_BYTES[] 
     {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //0-15 
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //16-31
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //32-47
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //48-63
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //64-79 (A,C,G)
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //80-95 (T)
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //96-111 (a,c,g)
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};//112-127 (t)

const char BASE_COMP_C[] 
     {'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //0-15  ga
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //16-31 rb
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //32-47 ag
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //48-63 e!
      'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N', //64-79 (A,C,G)
      'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N', //80-95 (T)
      'N','t','N','g','N','N','N','c','N','N','N','N','N','N','N','N', //96-111 (a,c,g)
      'N','N','N','N','a','N','N','N','N','N','N','N','N','N','N','N'};//112-127 (t)

const u8 BASE_COMP_B[] {3, 2, 1, 0};

template<typename T>
class ValArray: public std::valarray<T> {
    public:
    using std::valarray<T>::valarray;
    using super = std::valarray<T>;

    static constexpr T init_na() {
        return  std::numeric_limits<T>::has_quiet_NaN ? std::numeric_limits<T>::quiet_NaN() : std::numeric_limits<T>::min();
    }

    static constexpr T NA = init_na();//std::numeric_limits<T>::max();

    const T &at(size_t i) const {
        return super::operator[](i);
    }

    T &at(size_t i) {
        return super::operator[](i);
    }

    std::string to_string() const {
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
        } else if (super::size() > 0) {
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

    T mean() const {
        auto sum = super::sum();
        auto len = super::size();
        return super::sum() / super::size();
    }

    T stdv(T mean) const {
        auto deltas = (*this) - mean;
        return sqrt((deltas*deltas).sum() / super::size());
    }

    bool empty() const {
        return super::size() == 0;
    }

    T stdv() const {
        return stdv(mean());
    }

    static py::class_<ValArray> pybind(py::module_ &m, std::string suffix) {
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
        a.attr("NA") = py::cast(ValArray::NA);
        a.def("__getitem__", static_cast<T& (super::*)(size_t)> (&ValArray::at));
        //        py::vectorize(static_cast< KmerType (*) (const KmerTypePy &, u32)>(&Class::str_to_kmer)), 
        //a.def("__getitem__", &ValArray::at);
        a.def("__repr__", &ValArray::to_string);
        a.def("to_numpy", &ValArray::to_numpy);
        a.def("__len__", &super::size);
        a.def("__iter__", [](const ValArray &v) { 
                return py::make_iterator(std::begin(v), std::end(v)); 
            }, py::keep_alive<0, 1>());

        return a;
    }
};

//TODO turn into ArrayRef, could represent slices of any array?
template<typename T>
struct PyArray {

    py::buffer_info info;
    T *data;
    size_t size_;

    //PyArray(T *ptr, size_t length) :
    //    info { ptr, length },
    //    data { ptr },
    //    size_ { length } {
    //}

    PyArray() {}

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
    static py::class_<Subclass> pybind(py::module_ &m, std::string suffix) {
        auto name = "PyArray" + suffix;
        py::class_<Subclass> c(m, name.c_str(), py::buffer_protocol());
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
void pybind_array(py::module_ &m, std::string suffix) {
    ValArray<T>::pybind(m, suffix);
    PyArray<T>::pybind(m, suffix);
}

void pybind_arrays(py::module_ &m);

#endif
