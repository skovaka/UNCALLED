//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#ifndef __HDF5_TOOLS_HPP
#define __HDF5_TOOLS_HPP

#include <cassert>
#include <cstring>
#include <exception>
#include <functional>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <type_traits>

/// Original HDF5 C API.
#ifndef DOXY
#include <hdf5.h>
#endif

/// New C++ wrapper for the HDF5 C API.
namespace hdf5_tools
{

/// Exception class thrown by failed hdf5 operations.
class Exception
    : public std::exception
{
public:
    /// Ctor: from message @p msg.
    Exception(std::string const & msg) : _msg(active_path() + ": " + msg) {}
    /// Get exception message.
    char const * what() const noexcept { return _msg.c_str(); }
    /**
     * Active HDF5 operation path.
     * If set during a high-level HDF5 operation, it will be printed as part of any exception message.
     */
    static std::string & active_path()
    {
        static thread_local std::string _active_path;
        return _active_path;
    }
private:
    std::string _msg;
}; // class Exception

// Forward declaration
class Compound_Map;

/// Internal constructs.
namespace detail
{

/**
 * Compute offset of a struct member from a member pointer (runtime version).
 * @tparam T Struct type.
 * @tparam U Member type.
 */
template <typename T, typename U>
std::size_t offset_of(U T::* mem_ptr)
{
    return reinterpret_cast<std::size_t>(&(((T*)0)->*mem_ptr));
}

/**
 * TMF: Get HDF5 memory type for native numeric destination type.
 * Note HDF5 idiosyncrasy:
 *   Types such as H5T_NATIVE_INT are not constants(!?), so id() is not a constexpr.
 */
template <typename T>
struct get_mem_type
{
    /// Id
    static hid_t id()
    {
        return -1;
    }
};
#ifndef DOXY
// signed integral
template <> struct get_mem_type<char>               { static hid_t id() { return H5T_NATIVE_CHAR;    } };
template <> struct get_mem_type<short>              { static hid_t id() { return H5T_NATIVE_SHORT;   } };
template <> struct get_mem_type<int>                { static hid_t id() { return H5T_NATIVE_INT;     } };
template <> struct get_mem_type<long>               { static hid_t id() { return H5T_NATIVE_LONG;    } };
template <> struct get_mem_type<long long>          { static hid_t id() { return H5T_NATIVE_LLONG;   } };
// unsigned integral
template <> struct get_mem_type<unsigned char>      { static hid_t id() { return H5T_NATIVE_UCHAR;   } };
template <> struct get_mem_type<unsigned short>     { static hid_t id() { return H5T_NATIVE_USHORT;  } };
template <> struct get_mem_type<unsigned>           { static hid_t id() { return H5T_NATIVE_UINT;    } };
template <> struct get_mem_type<unsigned long>      { static hid_t id() { return H5T_NATIVE_ULONG;   } };
template <> struct get_mem_type<unsigned long long> { static hid_t id() { return H5T_NATIVE_ULLONG;  } };
// float
template <> struct get_mem_type<float>              { static hid_t id() { return H5T_NATIVE_FLOAT;   } };
template <> struct get_mem_type<double>             { static hid_t id() { return H5T_NATIVE_DOUBLE;  } };
template <> struct get_mem_type<long double>        { static hid_t id() { return H5T_NATIVE_LDOUBLE; } };
#endif

/**
 * TMF: Get class of memory type.
 * Results:
 * - 0: unknown
 * - 1: numeric (signed/unsigned integer or float)
 * - 2: fixed length string (char array)
 * - 3: variable length string (std::string)
 * - 4: struct/class
 */
template <typename T>
struct mem_type_class
{
#ifndef DOXY
    static int const value =
        std::conditional<std::is_integral<T>::value or std::is_floating_point<T>::value,
                         std::integral_constant<int, 1>,
                         typename std::conditional<std::is_class<T>::value,
                                                   std::integral_constant<int, 4>,
                                                   std::integral_constant<int, 0>>::type>::type::value;
#else
    /// Value
    static int const value = 0;
#endif
};
#ifndef DOXY
template <size_t Size>
struct mem_type_class<char[Size]>
{
    static int const value = 2;
};
template <size_t Size>
struct mem_type_class<char const[Size]>
{
    static int const value = 2;
};
template <size_t Size>
struct mem_type_class<std::array<char, Size>>
{
    static int const value = 2;
};
template <size_t Size>
struct mem_type_class<std::array<char const, Size>>
{
    static int const value = 2;
};
template <>
struct mem_type_class<std::string>
{
    static int const value = 3;
};
#endif

/**
 * HDF5 object holder.
 * Upon destruction, deallocate the held HDF5 object, check the HDF5 API return value, and throw exception on error.
 */
struct HDF_Object_Holder
{
    /// HDF5 object
    hid_t id;
    /// HDF5 object deallocator
    std::function<herr_t(hid_t)> dtor;

    /// Ctor: default
    HDF_Object_Holder()
        : id(0) {}
    /// Ctor: copy
    HDF_Object_Holder(HDF_Object_Holder const &) = delete;
    /// Ctor: move
    HDF_Object_Holder(HDF_Object_Holder && other)
        : id(0)
    {
        *this = std::move(other);
    }
    /// Ctor: from HDF5 object and destructor
    HDF_Object_Holder(hid_t _id, std::function<herr_t(hid_t)> _dtor)
        : id(_id), dtor(_dtor) {}
    /// Dtor: throw on HDF5 errors
    ~HDF_Object_Holder() noexcept(false)
    {
        if (id > 0)
        {
            if (dtor)
            {
                dtor(id);
            }
            id = 0;
        }
    }
    /// Asop: copy
    HDF_Object_Holder & operator = (HDF_Object_Holder const &) = delete;
    /// Asop: move
    HDF_Object_Holder & operator = (HDF_Object_Holder && other)
    {
        if (&other != this)
        {
            std::swap(id, other.id);
            std::swap(dtor, other.dtor);
        }
        return *this;
    }
}; // struct HDF_Object_Holder

/// Holder of static methods
struct Util
{
    /**
     * Create HDF5 string type.
     * @param sz If nonnegative, create fixed length string type of size @p sz;
     * else create variable length string type.
     */
    static HDF_Object_Holder
    make_str_type(long sz)
    {
        assert(sz != 0);
        HDF_Object_Holder res(
            wrap(H5Tcopy, H5T_C_S1),
            wrapped_closer(H5Tclose));
        size_t tmp = sz < 0? H5T_VARIABLE : sz;
        wrap(H5Tset_size, res.id, tmp);
        return res;
    } // make_str_type()

    /**
     * Wrapper of HDF5 calls.
     * Apply HDF5 call and check return value for validity.
     * @tparam Function HDF5 function type.
     * @tparam Args Variable arguments type.
     * @param f HDF5 function.
     * @param args Variable arguments.
     */
    template <typename Function, typename... Args>
    static typename std::result_of<Function(Args...)>::type
    wrap(Function && f, Args && ...args)
    {
        auto res = f(args...);
        auto const & f_info = get_fcn_info((void(*)())&f);
        if (not f_info.second((void*)&res)) throw Exception(std::string("error in ") + f_info.first);
        return res;
    }

    /**
     * Wrapped closer function.
     * Create and return a wrapped closer function.
     * @tparam Function HDF5 closer function type.
     * @param f HDF5 closer function.
     */
    template <typename Function>
    static std::function<herr_t(hid_t)>
    wrapped_closer(Function && f)
    {
        return [&] (hid_t id) { return wrap(f, id); };
    }

private:
    /**
     * Get name and return value checker for HDF5 function.
     * Maintain a static map of HDF5 functions and corresponding checks for return value validity.
     */
    static std::pair<char const *, std::function<bool(void *)>> const &
    get_fcn_info(void (*fcn_ptr)())
    {
        static std::map<void (*)(), std::pair<char const *, std::function<bool(void *)>>> const fcn_info_m =
            {
                { (void(*)())&H5Aclose,
                  { "H5Aclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Acreate2,
                  { "H5Acreate2",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Aexists_by_name,
                  { "H5Aexists_by_name",
                    [] (void * vp) { return *reinterpret_cast<htri_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Aget_name_by_idx,
                  { "H5Aget_name_by_idx",
                    [] (void * vp) { return *reinterpret_cast<ssize_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Aget_space,
                  { "H5Aget_space",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Aget_type,
                  { "H5Aget_type",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Aopen,
                  { "H5Aopen",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Aopen_by_name,
                  { "H5Aopen_by_name",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Aread,
                  { "H5Aread",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Awrite,
                  { "H5Awrite",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },

                { (void(*)())&H5Dclose,
                  { "H5Dclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Dcreate2,
                  { "H5Dcreate2",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Dget_space,
                  { "H5Dget_space",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Dget_type,
                  { "H5Dget_type",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Dopen,
                  { "H5Dopen",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Dread,
                  { "H5Dread",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Dvlen_reclaim,
                  { "H5Dvlen_reclaim",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Dwrite,
                  { "H5Dwrite",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },

                { (void(*)())&H5Gclose,
                  { "H5Gclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Gcreate2,
                  { "H5Gcreate2",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Gget_info,
                  { "H5Gget_info",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Gopen2,
                  { "H5Gopen2",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },

                { (void(*)())&H5Lexists,
                  { "H5Lexists",
                    [] (void * vp) { return *reinterpret_cast<htri_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Lget_name_by_idx,
                  { "H5Lget_name_by_idx",
                    [] (void * vp) { return *reinterpret_cast<ssize_t *>(vp) >= 0; }
                  }
                },

                { (void(*)())&H5Oclose,
                  { "H5Oclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Oexists_by_name,
                  { "H5Oexists_by_name",
                    [] (void * vp) { return *reinterpret_cast<htri_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Oget_info,
                  { "H5Oget_info",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Oopen,
                  { "H5Oopen",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },

                { (void(*)())&H5Pclose,
                  { "H5Pclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Pcreate,
                  { "H5Pcreate",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Pset_create_intermediate_group,
                  { "H5Pset_create_intermediate_group",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },

                { (void(*)())&H5Sclose,
                  { "H5Sclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Screate,
                  { "H5Screate",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Screate_simple,
                  { "H5Screate_simple",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Sget_simple_extent_dims,
                  { "H5Sget_simple_extent_dims",
                    [] (void * vp) { return *reinterpret_cast<int *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Sget_simple_extent_ndims,
                  { "H5Sget_simple_extent_ndims",
                    [] (void * vp) { return *reinterpret_cast<int *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Sget_simple_extent_type,
                  { "H5Sget_simple_extent_type",
                    [] (void * vp) { return *reinterpret_cast<H5S_class_t *>(vp) != H5S_NO_CLASS; }
                  }
                },

                { (void(*)())&H5Tclose,
                  { "H5Tclose",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Tcopy,
                  { "H5Tcopy",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Tcreate,
                  { "H5Tcreate",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Tget_class,
                  { "H5Tget_class",
                    [] (void * vp) { return *reinterpret_cast<H5T_class_t *>(vp) != H5T_NO_CLASS; }
                  }
                },
                { (void(*)())&H5Tget_cset,
                  { "H5Tget_cset",
                    [] (void * vp) { return *reinterpret_cast<H5T_cset_t *>(vp) != H5T_CSET_ERROR; }
                  }
                },
                { (void(*)())&H5Tget_member_index,
                  { "H5Tget_member_index",
                    [] (void * vp) { return *reinterpret_cast<int *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Tget_member_name,
                  { "H5Tget_member_name",
                    [] (void * vp) { return *reinterpret_cast<char* *>(vp) != nullptr; }
                  }
                },
                { (void(*)())&H5Tget_member_type,
                  { "H5Tget_member_type",
                    [] (void * vp) { return *reinterpret_cast<hid_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Tget_nmembers,
                  { "H5Tget_nmembers",
                    [] (void * vp) { return *reinterpret_cast<int *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Tget_sign,
                  { "H5Tget_sign",
                    [] (void * vp) { return *reinterpret_cast<H5T_sign_t *>(vp) != H5T_SGN_ERROR; }
                  }
                },
                { (void(*)())&H5Tget_size,
                  { "H5Tget_size",
                    [] (void * vp) { return *reinterpret_cast<size_t *>(vp) > 0; }
                  }
                },
                { (void(*)())&H5Tinsert,
                  { "H5Tinsert",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Tis_variable_str,
                  { "H5Tis_variable_str",
                    [] (void * vp) { return *reinterpret_cast<htri_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Tset_cset,
                  { "H5Tset_cset",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
                { (void(*)())&H5Tset_size,
                  { "H5Tset_size",
                    [] (void * vp) { return *reinterpret_cast<herr_t *>(vp) >= 0; }
                  }
                },
            };
        return fcn_info_m.at(fcn_ptr);
    }

}; // struct Util

/**
 * Description of a compound member.
 * Only for numeric, string, and struct member types.
 */
struct Compound_Member_Description
{
public:
    /**
     * Ctor: numeric compound member
     * @param _name Compound member name
     * @param _offset Compound member offset
     * @param _numeric_type_id HDF5 numeric type
     */
    Compound_Member_Description(std::string const & _name, size_t _offset, hid_t _numeric_type_id)
        : type(numeric),
          name(_name),
          offset(_offset),
          numeric_type_id(_numeric_type_id) {}
    /**
     * Ctor: fixed-length string compound member
     * @param _name Compound member name
     * @param _offset Compound member offset
     * @param _char_array_size Fixed length string size
     */
    Compound_Member_Description(std::string const & _name, size_t _offset, size_t _char_array_size)
        : type(char_array),
          name(_name),
          offset(_offset),
          char_array_size(_char_array_size) {}
    /**
     * Ctor: variable-length string compound member
     * @param _name Compound member name
     * @param _offset Compound member offset
     */
    Compound_Member_Description(std::string const & _name, size_t _offset)
        : type(string),
          name(_name),
          offset(_offset) {}
    /**
     * Ctor: struct compound member
     * @param _name Compound member name
     * @param _offset Compound member offset
     * @param _compound_map_ptr Inner compound map
     * @param _compound_size Inner compound size
     */
    Compound_Member_Description(std::string const & _name, size_t _offset,
                                Compound_Map const * _compound_map_ptr, size_t _compound_size)
        : type(compound),
          name(_name),
          offset(_offset),
          compound_map_ptr(_compound_map_ptr),
          compound_size(_compound_size) {}

    /// Is member a numeric
    bool is_numeric()    const { return type == numeric;    }
    /// Is member a fixed-length string
    bool is_char_array() const { return type == char_array; }
    /// Is member a variable-length string
    bool is_string()     const { return type == string;     }
    /// Is member a compound
    bool is_compound()   const { return type == compound;   }

    /// Create and return a holder for the HDF5 type of this compound member.
    HDF_Object_Holder get_type() const
    {
        assert(not is_compound());
        HDF_Object_Holder res;
        if (is_numeric())
        {
            res = HDF_Object_Holder(numeric_type_id, nullptr);
        }
        else if (is_char_array())
        {
            res = Util::make_str_type(char_array_size);
        }
        else if (is_string())
        {
            res = Util::make_str_type(-1);
        }
        return res;
    } // get_type()

    /// Write compound member description to stream
    friend std::ostream & operator << (std::ostream & os, Compound_Member_Description const & e)
    {
        os << "(&=" << (void*)&e
           << ",name=\"" << e.name
           << "\",type=" << (e.is_numeric()
                             ? "numeric"
                             : (e.is_char_array()
                                ? "char_array"
                                : (e.is_string()
                                   ? "string" : "compound")))
           << ",offset=" << e.offset << ")";
        return os;
    }

    /// Compound member type enum
    enum member_type
    {
        numeric,
        char_array,
        string,
        compound
    };

    /// Compound member type
    member_type type;
    /// Compound member name
    std::string name;
    /// Compound member offset
    size_t offset;
    union
    {
        /// For numeric types: HDF5 type
        hid_t numeric_type_id;
        /// For string types: size of string type, 0 if variable-length.
        size_t char_array_size;
        /// For compound types: inner compound map
        Compound_Map const * compound_map_ptr;
    };
    /// For compound types: inner compound size
    size_t compound_size;
}; // struct Compound_Member_Description

} // namespace detail

/// A map of struct fields to tags that is used to read compound datatypes
class Compound_Map
{
public:
    /// Ctor: default
    Compound_Map() = default;
    /// Ctor: copy
    Compound_Map(Compound_Map const &) = delete;
    /// Ctor: move
    Compound_Map(Compound_Map &&) = default;
    /// Asop: copy
    Compound_Map & operator = (Compound_Map const &) = delete;
    /// Asop: move
    Compound_Map & operator = (Compound_Map &&) = default;

    /**
     * Add numeric/string compound member.
     * @tparam T struct type
     * @tparam U member type
     * @param name HDF5 compound member name
     * @param mem_ptr struct member pointer
     */
    template <typename T, typename U>
    void add_member(std::string const & name, U T::* mem_ptr)
    {
        static_assert(detail::mem_type_class<U>::value == 1
                      or detail::mem_type_class<U>::value == 2
                      or detail::mem_type_class<U>::value == 3,
                      "add_member(name, mem_ptr) overload expects numerical or string types only");
        if (detail::mem_type_class<U>::value == 1)
        {
            _members.emplace_back(name, detail::offset_of(mem_ptr), detail::get_mem_type<U>::id());
        }
        else if (detail::mem_type_class<U>::value == 2)
        {
            _members.emplace_back(name, detail::offset_of(mem_ptr), sizeof(U));
        }
        else if (detail::mem_type_class<U>::value == 3)
        {
            _members.emplace_back(name, detail::offset_of(mem_ptr));
        }
    }

    /**
     * Add compound compound member.
     * @tparam T struct type
     * @tparam U member type
     * @param name HDF5 compound member name
     * @param mem_ptr struct member pointer
     * @param compound_map_ptr inner compound map
     */
    template <typename T, typename U>
    void add_member(std::string const & name, U T::* mem_ptr, Compound_Map const * compound_map_ptr)
    {
        static_assert(detail::mem_type_class<U>::value == 4,
                      "add_member(name, mem_ptr, compound_map_ptr) overload expects class types only");
        _members.emplace_back(name, detail::offset_of(mem_ptr), compound_map_ptr, sizeof(U));
    }

    /// Get list of compound members.
    std::vector<detail::Compound_Member_Description> const & members() const { return _members; }

    /// Type returned by @p get_member_ptr_list().
    typedef std::deque<std::pair<std::deque<detail::Compound_Member_Description const *>,
                                 size_t>> member_ptr_list_type;
    /**
     * Get list of inner non-compound member types.
     * If members of the top-level compound are also compounds, they are traversed recursively.
     * @return A list of pairs; first: list of member ptrs followed; second: absolute offset
     */
    member_ptr_list_type get_member_ptr_list() const
    {
        member_ptr_list_type res;
        for (auto const & e : members())
        {
            member_ptr_list_type::value_type p;
            if (not e.is_compound())
            {
                member_ptr_list_type::value_type p;
                p.first = { &e };
                p.second = e.offset;
                res.emplace_back(std::move(p));
            }
            else
            {
                auto tmp = e.compound_map_ptr->get_member_ptr_list();
                for (auto & tmp_p : tmp)
                {
                    member_ptr_list_type::value_type p;
                    p.first = std::move(tmp_p.first);
                    p.first.push_front(&e);
                    p.second = tmp_p.second + e.offset;
                    res.emplace_back(std::move(p));
                }
            }
        }
        return res;
    }

    /**
     * Create an HDF5 compound type for this map.
     * @param compound_size Extrenally-tracked compound size
     * @param selector If empty, use all elements; if not empty, use only elements that pass selection
     * @param fill If true, type offsets follow compound map offsets, allowing for gaps;
     * if false: type offsets are minimal values required to fit members
     */
    detail::HDF_Object_Holder build_type(
        size_t compound_size,
        std::function<bool(detail::Compound_Member_Description const &)> selector = nullptr,
        bool fill = true) const
    {
        //std::clog << "===== build_type (" << (void*)this << ") start" << std::endl;
        std::deque<std::tuple<std::string, detail::HDF_Object_Holder, size_t>> stype_id_holder_l;
        size_t compressed_size = 0;
        for (auto const & e : members())
        {
            detail::HDF_Object_Holder stype_id_holder;
            if (selector and not e.is_compound() and not selector(e)) continue;
            if (not e.is_compound())
            {
                stype_id_holder = e.get_type();
            }
            else
            {
                stype_id_holder = e.compound_map_ptr->build_type(e.compound_size, selector, fill);
            }
            if (stype_id_holder.id > 0)
            {
                stype_id_holder_l.emplace_back(
                    std::string(e.name),
                    std::move(stype_id_holder),
                    fill? e.offset : compressed_size);
                compressed_size += H5Tget_size(std::get<1>(stype_id_holder_l.back()).id);
            }
        }
        if (stype_id_holder_l.empty())
        {
            //std::clog << "===== build_type (" << (void*)this << ") empty" << std::endl;
            return detail::HDF_Object_Holder();
        }
        //std::clog << "===== build_type (" << (void*)this << ") compound size: " << (fill? compound_size : compressed_size) << std::endl;
        detail::HDF_Object_Holder res(
            detail::Util::wrap(H5Tcreate, H5T_COMPOUND, fill? compound_size : compressed_size),
            detail::Util::wrapped_closer(H5Tclose));
        for (auto const & t : stype_id_holder_l)
        {
            //std::clog << "===== build_type (" << (void*)this << ") adding name=\"" << std::get<0>(t) << "\", offset=" << std::get<2>(t) << std::endl;
            detail::Util::wrap(H5Tinsert, res.id, std::get<0>(t).c_str(), std::get<2>(t), std::get<1>(t).id);
        }
        //std::clog << "===== build_type (" << (void*)this << ") end" << std::endl;
        return res;
    }

    /**
     * Create a flat HDF5 compound type for members in the given list
     * @param l List of compound member descriptions
     * @param id If non-zero, HDF5 id of compound type to extend
     */
    static detail::HDF_Object_Holder build_flat_type(
        member_ptr_list_type::value_type::first_type const & l, hid_t id = 0)
    {
        detail::HDF_Object_Holder res;
        size_t sz = 0;
        for (auto it = l.rbegin(); it != l.rend(); ++it)
        {
            auto const & e = **it;
            assert((it == l.rbegin()) == (not e.is_compound()));
            assert((it == l.rbegin()) == (res.id == 0));
            assert((it == l.rbegin()) == (sz == 0));
            if (it == l.rbegin())
            {
                if (id == 0)
                {
                    res = e.get_type();
                }
                else
                {
                    res = detail::HDF_Object_Holder(
                        detail::Util::wrap(H5Tcopy, id),
                        detail::Util::wrapped_closer(H5Tclose));
                }
                sz = detail::Util::wrap(H5Tget_size, res.id);
            }
            detail::HDF_Object_Holder tmp(
                detail::Util::wrap(H5Tcreate, H5T_COMPOUND, sz),
                detail::Util::wrapped_closer(H5Tclose));
            detail::Util::wrap(H5Tinsert, tmp.id, e.name.c_str(), 0, res.id);
            std::swap(res, tmp);
        }
        return res;
    }

    /**
     * Extract inner compound member type.
     * Given an existing compound type and a list of member compound pointers to follow,
     * extract inner compound type.
     * @param id HDF5 compound type
     * @param l List of inner compound members to follow (all but the last must be compounds)
     */
    static detail::HDF_Object_Holder get_compound_member(
        hid_t id, member_ptr_list_type::value_type::first_type const & l)
    {
        detail::HDF_Object_Holder res(
            detail::Util::wrap(H5Tcopy, id),
            detail::Util::wrapped_closer(H5Tclose));
        for (auto it = l.begin(); it != l.end(); ++it)
        {
            auto const & e = **it;
            assert(detail::Util::wrap(H5Tget_class, res.id) == H5T_COMPOUND);
            unsigned idx = detail::Util::wrap(H5Tget_member_index, res.id, e.name.c_str());
            detail::HDF_Object_Holder tmp(
                detail::Util::wrap(H5Tget_member_type, res.id, idx),
                detail::Util::wrapped_closer(H5Tclose));
            std::swap(res, tmp);
        }
        assert(detail::Util::wrap(H5Tget_class, res.id) != H5T_COMPOUND);
        return res;
    }

private:
    std::vector<detail::Compound_Member_Description> _members;
}; // Compound_Map

namespace detail
{

/**
 * Manage HDF5 file object access during a read operation.
 * Opens a file object, and retrieves various attributes.
 * Note: called by @p Reader; level 2 below @p File::read().
 */
struct Reader_Base
{
    /**
     * Ctor: from HDF5 group and name.
     * @param grp_id HDF5 group id
     * @param name HDF5 object name (dataset or attribute)
     */
    Reader_Base(hid_t grp_id, std::string const & name)
    {
        int status = Util::wrap(H5Aexists_by_name, grp_id, ".", name.c_str(), H5P_DEFAULT);
        is_ds = status == 0;
        if (is_ds)
        {
            obj_id_holder = HDF_Object_Holder(
                Util::wrap(H5Dopen, grp_id, name.c_str(), H5P_DEFAULT),
                Util::wrapped_closer(H5Dclose));
            dspace_id_holder = HDF_Object_Holder(
                Util::wrap(H5Dget_space, obj_id_holder.id),
                Util::wrapped_closer(H5Sclose));
            file_dtype_id_holder = HDF_Object_Holder(
                Util::wrap(H5Dget_type, obj_id_holder.id),
                Util::wrapped_closer(H5Tclose));
            reader = [&] (hid_t mem_dtype_id, void * dest) {
                return Util::wrap(H5Dread, obj_id_holder.id, mem_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dest);
            };
        }
        else
        {
            obj_id_holder = HDF_Object_Holder(
                Util::wrap(H5Aopen, grp_id, name.c_str(), H5P_DEFAULT),
                Util::wrapped_closer(H5Aclose));
            dspace_id_holder = HDF_Object_Holder(
                Util::wrap(H5Aget_space, obj_id_holder.id),
                Util::wrapped_closer(H5Sclose));
            file_dtype_id_holder = HDF_Object_Holder(
                Util::wrap(H5Aget_type, obj_id_holder.id),
                Util::wrapped_closer(H5Tclose));
            reader = [&] (hid_t mem_dtype_id, void * dest) {
                return Util::wrap(H5Aread, obj_id_holder.id, mem_dtype_id, dest);
            };
        }
        // dataspace class and size
        dspace_class = Util::wrap(H5Sget_simple_extent_type, dspace_id_holder.id);
        if (dspace_class == H5S_SCALAR)
        {
            dspace_size = 1;
        }
        else if (dspace_class == H5S_SIMPLE)
        {
            auto ndims = Util::wrap(H5Sget_simple_extent_ndims, dspace_id_holder.id);
            if (ndims != 1) throw Exception("reading multi-dimensional extents is not supported");
            hsize_t tmp;
            Util::wrap(H5Sget_simple_extent_dims, dspace_id_holder.id, &tmp, nullptr);
            dspace_size = tmp;
        }
        else
        {
            throw Exception("reading dataspaces other than SCALAR and SIMPLE is not supported");
        }
        // datatype class
        file_dtype_class = Util::wrap(H5Tget_class, file_dtype_id_holder.id);
        if (file_dtype_class == H5T_STRING)
        {
            file_dtype_is_vlen_str = Util::wrap(H5Tis_variable_str, file_dtype_id_holder.id);
        }
        else
        {
            file_dtype_is_vlen_str = false;
        }
        // datatype size
        file_dtype_size = Util::wrap(H5Tget_size, file_dtype_id_holder.id);
    }
    /// Object holder
    HDF_Object_Holder obj_id_holder;
    /// Dataspace holder
    HDF_Object_Holder dspace_id_holder;
    /// Datatype holder
    HDF_Object_Holder file_dtype_id_holder;
    /// Reader function
    std::function<void(hid_t, void *)> reader;
    /// Dataspace class
    H5S_class_t dspace_class;
    /// Dataspace size
    size_t dspace_size;
    /// Datatype class
    H5T_class_t file_dtype_class;
    /// Variable-length string flag
    htri_t file_dtype_is_vlen_str;
    /// Datatype size
    size_t file_dtype_size;
    /// Is-dataset flag
    bool is_ds;
}; // struct Reader_Base

/**
 * Read a vector of strings.
 */
struct String_Reader
{
    /**
     * Perform read.
     * Given a file object and an optional inner member pointer selector,
     * read extent into vector of strings.
     * The function converts numeric and fixed/variable-length string HDF5 types to string.
     * @param reader_base File object read manager.
     * @param mptr_l_ptr Optional inner compound member selector. If file object is a compound,
     * the list should specify which inner compound members to follow to get to a base
     * non-compound member to read.
     */
    std::vector<std::string> operator () (
        Reader_Base const & reader_base,
        Compound_Map::member_ptr_list_type::value_type::first_type const * mptr_l_ptr = nullptr) const
    {
        std::vector<std::string> res(reader_base.dspace_size);
        assert((mptr_l_ptr != nullptr) == (reader_base.file_dtype_class == H5T_COMPOUND));
        HDF_Object_Holder file_stype_id_holder;
        hid_t file_stype_id = 0;
        if (reader_base.file_dtype_class == H5T_COMPOUND)
        {
            file_stype_id_holder = Compound_Map::get_compound_member(
                reader_base.file_dtype_id_holder.id,
                *mptr_l_ptr);
            file_stype_id = file_stype_id_holder.id;
        }
        else
        {
            file_stype_id = reader_base.file_dtype_id_holder.id;
        }
        auto mem_type_wrapper = [&] (HDF_Object_Holder && id_holder) {
            HDF_Object_Holder tmp(std::move(id_holder));
            return (mptr_l_ptr != nullptr
                    ? Compound_Map::build_flat_type(*mptr_l_ptr, tmp.id)
                    : std::move(tmp));
        };
        assert(Util::wrap(H5Tget_class, file_stype_id) != H5T_COMPOUND);
        auto file_stype_class = Util::wrap(H5Tget_class, file_stype_id);
        HDF_Object_Holder mem_dtype_id_holder;
        if (file_stype_class == H5T_STRING) // stored as a string
        {
            auto file_stype_cset = Util::wrap(H5Tget_cset, file_stype_id);
            if (Util::wrap(H5Tis_variable_str, file_stype_id)) // stored as a varlen string
            {
                // compute mem_type
                auto mem_stype_id_holder = Util::make_str_type(-1);
                Util::wrap(H5Tset_cset, mem_stype_id_holder.id, file_stype_cset);
                mem_dtype_id_holder = mem_type_wrapper(std::move(mem_stype_id_holder));
                // prepare buffer to receive data
                std::vector<char *> charptr_buff(res.size(), nullptr);
                // perform the read
                reader_base.reader(mem_dtype_id_holder.id, charptr_buff.data());
                // transfer strings to destination
                for (size_t i = 0; i < res.size(); ++i)
                {
                    if (not charptr_buff[i]) throw Exception("read did not fill buffer");
                    res[i] = charptr_buff[i];
                }
                // reclaim memory allocated by libhdf5
                Util::wrap(H5Dvlen_reclaim, mem_dtype_id_holder.id, reader_base.dspace_id_holder.id,
                           H5P_DEFAULT, charptr_buff.data());
            }
            else // stored as a fixlen string
            {
                // compute mem_type
                size_t file_stype_size = Util::wrap(H5Tget_size, file_stype_id);
                auto mem_stype_id_holder = Util::make_str_type(file_stype_size + 1);
                Util::wrap(H5Tset_cset, mem_stype_id_holder.id, file_stype_cset);
                mem_dtype_id_holder = mem_type_wrapper(std::move(mem_stype_id_holder));
                // prepare buffer to receieve data
                std::vector<char> char_buff(res.size() * (file_stype_size + 1), '\0');
                // perform the read
                reader_base.reader(mem_dtype_id_holder.id, char_buff.data());
                // transfer strings to destination
                for (size_t i = 0; i < res.size(); ++i)
                {
                    res[i] = std::string(&char_buff[i * (file_stype_size + 1)], file_stype_size);
                    // trim trailing '\0'-s
                    while (not res[i].empty() and res[i].back() == '\0')
                    {
                        res[i].resize(res[i].size() - 1);
                    }
                }
            }
        }
        else if (file_stype_class == H5T_INTEGER) // stored as an integer
        {
            if (Util::wrap(H5Tget_sign, file_stype_id) == H5T_SGN_NONE) // stored as an unsigned integer
            {
                // compute mem_type
                mem_dtype_id_holder = mem_type_wrapper(
                    HDF_Object_Holder(get_mem_type<unsigned long long>::id(), nullptr));
                // prepare buffer to read data
                std::vector<unsigned long long> ull_buff(res.size());
                // perform the read
                reader_base.reader(mem_dtype_id_holder.id, ull_buff.data());
                // transfer to destination
                for (size_t i = 0; i < res.size(); ++i)
                {
                    std::ostringstream oss;
                    oss << ull_buff[i];
                    res[i] = oss.str();
                }
            }
            else // stored as a signed integer
            {
                // compute mem_type
                mem_dtype_id_holder = mem_type_wrapper(
                    HDF_Object_Holder(get_mem_type<long long>::id(), nullptr));
                // prepare buffer to read data
                std::vector<long long> ll_buff(res.size());
                // perform the read
                reader_base.reader(mem_dtype_id_holder.id, ll_buff.data());
                // transfer to destination
                for (size_t i = 0; i < res.size(); ++i)
                {
                    std::ostringstream oss;
                    oss << ll_buff[i];
                    res[i] = oss.str();
                }
            }
        }
        else if (file_stype_class == H5T_FLOAT) // stored as a float
        {
            // compute mem_type
            mem_dtype_id_holder = mem_type_wrapper(
                HDF_Object_Holder(get_mem_type<double>::id(), nullptr));
            // prepare buffer to read data
            std::vector<double> d_buff(res.size());
            // perform the read
            reader_base.reader(mem_dtype_id_holder.id, d_buff.data());
            // transfer to destination
            for (size_t i = 0; i < res.size(); ++i)
            {
                std::ostringstream oss;
                oss << d_buff[i];
                res[i] = oss.str();
            }
        }
        return res;
    }
}; // struct String_Reader

/**
 * Construct memory data type and invoke @p Reader_Base reader.
 * Note: called by @p Reader; level 2 below @p File::read().
 */
template <int Mem_Type_Class, typename Data_Type>
struct Reader_Helper;

/// @details @em Specialization_Numeric: @p Mem_Type_Class = 1.
#ifndef DOXY
template <typename Data_Type>
struct Reader_Helper<1, Data_Type>
#else
template <int Mem_Type_Class, typename Data_Type>
struct Reader_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Numeric.
     * Obtain memory type from @p get_mem_type.
     * @param reader_base File object read manager
     * @param out Destination
     */
    void operator () (Reader_Base const & reader_base, Data_Type * out) const
    {
        assert(std::is_integral<Data_Type>::value or std::is_floating_point<Data_Type>::value);
        hid_t mem_dtype_id = get_mem_type<Data_Type>::id();
        reader_base.reader(mem_dtype_id, out);
    }
};

/// @details @em Specialization_Fix_Len_String: @p Mem_Type_Class = 2.
#ifndef DOXY
template <typename Data_Type>
struct Reader_Helper<2, Data_Type>
#else
template <int Mem_Type_Class, typename Data_Type>
struct Reader_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Fix_Len_String.
     * If the file datatype is also fixed-length, use @p Reader_Base reader directly.
     * Else, use @p String_Reader.
     * @param reader_base File object read manager
     * @param out Destination
     */
    void operator () (Reader_Base const & reader_base, Data_Type * out) const
    {
        if (reader_base.file_dtype_class == H5T_STRING
            and not reader_base.file_dtype_is_vlen_str)
        {
            HDF_Object_Holder mem_dtype_id_holder(Util::make_str_type(sizeof(Data_Type)));
            auto file_dtype_cset = Util::wrap(H5Tget_cset, reader_base.file_dtype_id_holder.id);
            Util::wrap(H5Tset_cset, mem_dtype_id_holder.id, file_dtype_cset);
            reader_base.reader(mem_dtype_id_holder.id, out);
        }
        else // conversion needed
        {
            auto tmp = String_Reader()(reader_base);
            for (size_t i = 0; i < tmp.size(); ++i)
            {
                std::memset(&out[i][0], '\0', sizeof(Data_Type));
                std::memcpy(&out[i][0], tmp[i].data(), std::min(tmp[i].size(), sizeof(Data_Type) - 1));
            }
        }
    }
};

/// @details @em Specialization_Var_Len_String: @p Mem_Type_Class = 3.
#ifndef DOXY
template <typename Data_Type>
struct Reader_Helper<3, Data_Type>
#else
template <int Mem_Type_Class, typename Data_Type>
struct Reader_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Var_Len_String.
     * Use @p String_Reader.
     * @param reader_base File object read manager
     * @param out Destination
     */
    void operator () (Reader_Base const & reader_base, Data_Type * out) const
    {
        static_assert(std::is_same<Data_Type, std::string>::value, "Data_Type not std::string");
        auto tmp = String_Reader()(reader_base);
        for (size_t i = 0; i < tmp.size(); ++i)
        {
            std::swap(out[i], tmp[i]);
        }
    }
};

/// @details @em Specialization_Compound: @p Mem_Type_Class = 4.
#ifndef DOXY
template <typename Data_Type>
struct Reader_Helper<4, Data_Type>
#else
template <int Mem_Type_Class, typename Data_Type>
struct Reader_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Compound.
     * Go through the list of compound members;
     * read all of the ones that don't need conversion in one go;
     * read the remaining ones one at a time.
     * @param reader_base File object read manager
     * @param out Destination
     * @param cm Compound map
     */
    void operator () (Reader_Base const & reader_base, Data_Type * out, Compound_Map const & cm) const
    {
        // get member list
        auto mptr_l = cm.get_member_ptr_list();
        // go through members, check they exist, decide if they need conversion
        std::set<detail::Compound_Member_Description const *> conversion_needed_s;
        for (auto const & p : mptr_l)
        {
            HDF_Object_Holder file_stype_id_holder(
                Compound_Map::get_compound_member(reader_base.file_dtype_id_holder.id, p.first));
            if (p.first.back()->is_string()
                or (p.first.back()->is_char_array()
                    and Util::wrap(H5Tget_class, file_stype_id_holder.id) == H5T_STRING
                    and Util::wrap(H5Tis_variable_str, file_stype_id_holder.id)))
            {
                conversion_needed_s.insert(p.first.back());
            }
        }
        // read all members that do not need conversion all-at-once
        auto implicit_conversion = [&] (detail::Compound_Member_Description const & e) {
            return conversion_needed_s.count(&e) == 0;
        };
        HDF_Object_Holder mem_dtype_id_holder(cm.build_type(sizeof(Data_Type), implicit_conversion, true));
        if (mem_dtype_id_holder.id > 0)
        {
            reader_base.reader(mem_dtype_id_holder.id, out);
        }
        // read members that need conversion one-by-one
        for (auto const & p : mptr_l)
        {
            auto const & e = *p.first.back();
            if (implicit_conversion(e)) continue;
            // read member into vector of strings
            auto tmp = String_Reader()(reader_base, &p.first);
            assert(tmp.size() == reader_base.dspace_size);
            // write it to destination
            assert(e.is_char_array() or e.is_string());
            if (e.is_char_array())
            {
                for (size_t i = 0; i < tmp.size(); ++i)
                {
                    std::memset(reinterpret_cast<char *>(&out[i]) + p.second, '\0', e.char_array_size);
                    std::memcpy(reinterpret_cast<char *>(&out[i]) + p.second,
                                tmp[i].data(),
                                std::min(tmp[i].size(), e.char_array_size - 1));
                }
            }
            else if (e.is_string())
            {
                for (size_t i = 0; i < tmp.size(); ++i)
                {
                    std::swap(
                        *reinterpret_cast<std::string *>(reinterpret_cast<char *>(&out[i]) + p.second),
                        tmp[i]);
                }
            }
        }
    }
};

/**
 * Prepare destination for a read operation.
 * Branch statically on the destination data type @p Dest_Type.
 * Use @p Reader_Base to open file object, then pass control to @p Reader_Helper.
 * Note: called by @p File::read.
 */
template <typename Dest_Type>
struct Reader;

/**
 * @struct Reader
 * @details @em Specialization_Default: @p Dest_Type = @p Data_Type.
 */
#ifndef DOXY
template <typename Data_Type>
struct Reader
#else
template <typename Dest_Type>
struct Reader
#endif
{
    /**
     * Functor operator.
     * @em Specialization_Default.
     * Check file object is scalar (dspace_size == 1).
     * Quirk: If reading a @p std::string, allow for file object to be an extent
     * of size 1 characters.
     * @param grp_id HDF5 group object
     * @param name HDF5 object name to read
     * @param out Destination (single address)
     * @param args Optional reading arguments passed to @p Reader_Helper
     */
    template <typename ...Args>
    void operator () (hid_t grp_id, std::string const & name,
                      Data_Type & out,
                      Args && ...args) const
    {
        Reader_Base reader_base(grp_id, name);
        if (reader_base.dspace_size == 1)
        {
            Reader_Helper<mem_type_class<Data_Type>::value, Data_Type>()(
                reader_base, &out, std::forward<Args>(args)...);
        }
        else if (std::is_same<Data_Type, std::string>::value
                 and reader_base.file_dtype_class == H5T_STRING
                 and not reader_base.file_dtype_is_vlen_str
                 and reader_base.file_dtype_size == 1)
        {
            std::vector<std::array<char, 1>> char_buff(reader_base.dspace_size);
            Reader_Helper<2, std::array<char, 1>>()(
                reader_base, char_buff.data(), std::forward<Args>(args)...);
            reinterpret_cast<std::string &>(out).assign(&char_buff[0][0], reader_base.dspace_size);
        }
        else
        {
            throw Exception("reading scalar, but dataspace size is not 1");
        }
    }
};
/**
 * @struct Reader
 * @details @em Specialization_Vector: @p Dest_Type = @p std::vector<Data_Type>.
 */
#ifndef DOXY
template <typename Data_Type>
struct Reader<std::vector<Data_Type>>
#else
template <typename Dest_Type>
struct Reader
#endif
{
    /**
     * Functor operator.
     * @em Specialization_Vector.
     * Clear and resize destination vector.
     * @param grp_id HDF5 group object
     * @param name HDF5 object name to read
     * @param out Destination (vector)
     * @param args Optional reading arguments passed to @p Reader_Helper
     */
    template <typename ...Args>
    void operator () (hid_t grp_id, std::string const & name,
                      std::vector<Data_Type> & out,
                      Args && ...args) const
    {
        Reader_Base reader_base(grp_id, name);
        out.clear();
        out.resize(reader_base.dspace_size);
        Reader_Helper<mem_type_class<Data_Type>::value, Data_Type>()(
            reader_base, out.data(), std::forward<Args>(args)...);
    }
};

/**
 * Holder of HDF5 access methods used during a write operation.
 */
struct Writer_Base
{
    /**
     * Create a dataset/attribute.
     * @param grp_id Parent group object.
     * @param loc_name dataset/attribute name.
     * @param as_ds Flag; create dataset iff true.
     * @param dspace_id HDF5 dataspace.
     * @param file_dtype_id HDF5 file datatype.
     */
    static HDF_Object_Holder create(hid_t grp_id, std::string const & loc_name, bool as_ds,
                                    hid_t dspace_id, hid_t file_dtype_id)
    {
        HDF_Object_Holder obj_id_holder;
        if (as_ds)
        {
            obj_id_holder = HDF_Object_Holder(
                Util::wrap(H5Dcreate2, grp_id, loc_name.c_str(), file_dtype_id, dspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                Util::wrapped_closer(H5Dclose));
        }
        else
        {
            obj_id_holder = HDF_Object_Holder(
                Util::wrap(H5Acreate2, grp_id, loc_name.c_str(), file_dtype_id, dspace_id,
                           H5P_DEFAULT, H5P_DEFAULT),
                Util::wrapped_closer(H5Aclose));
        }
        return obj_id_holder;
    }
    /**
     * Write dataset/attribute to file.
     * @param obj_id Destination object.
     * @param as_ds Flag; write dataset iff true.
     * @param mem_dtype_id HDF5 memory datatype.
     * @param in Source memory address.
     */
    static void write(hid_t obj_id, bool as_ds, hid_t mem_dtype_id, void const * in)
    {
        if (as_ds)
        {
            Util::wrap(H5Dwrite, obj_id, mem_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, in);
        }
        else
        {
            Util::wrap(H5Awrite, obj_id, mem_dtype_id, in);
        }
    }
    /**
     * Create and write dataset/attribute.
     * @param grp_id Parent group object.
     * @param loc_name dataset/attribute name.
     * @param as_ds Flag; create dataset iff true.
     * @param dspace_id HDF5 dataspace.
     * @param mem_dtype_id HDF5 memory datatype.
     * @param file_dtype_id HDF5 file datatype.
     * @param in Source memory address.
     */
    static void create_and_write(hid_t grp_id, std::string const & loc_name, bool as_ds,
                                 hid_t dspace_id, hid_t mem_dtype_id, hid_t file_dtype_id,
                                 void const * in)
    {
        HDF_Object_Holder obj_id_holder(create(grp_id, loc_name, as_ds, dspace_id, file_dtype_id));
        write(obj_id_holder.id, as_ds, mem_dtype_id, in);
    }
}; // struct Writer_Base


/**
 * Manage calls to Write_Base methods.
 * Note: called by @p Writer; level 2 below @p File::write().
 */
template <int, typename>
struct Writer_Helper;

/// @details @em Specialization_Numeric: @p Mem_Type_Class = 1.
#ifndef DOXY
template <typename In_Data_Type>
struct Writer_Helper<1, In_Data_Type>
#else
template <int Mem_Type_Class, typename In_Data_Type>
struct Writer_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Numeric.
     * Obtain memory type from @p get_mem_type.
     * @param grp_id Parent group object.
     * @param loc_name dataset/attribute name.
     * @param as_ds Flag; create dataset iff true.
     * @param dspace_id HDF5 dataspace.
     * @param in Source memory address.
     * @param file_dtype_id HDF5 file datatype;
     * if 0, use @p mem_dtype_id.
     */
    void operator () (hid_t grp_id, std::string const & loc_name,
                      bool as_ds, hid_t dspace_id, size_t,
                      In_Data_Type const * in, hid_t file_dtype_id = 0) const
    {
        assert(std::is_integral<In_Data_Type>::value or std::is_floating_point<In_Data_Type>::value);
        hid_t mem_dtype_id = get_mem_type<In_Data_Type>::id();
        if (file_dtype_id == 0)
        {
            file_dtype_id = mem_dtype_id;
        }
        Writer_Base::create_and_write(
            grp_id, loc_name, as_ds,
            dspace_id, mem_dtype_id, file_dtype_id,
            in);
    }
}; // struct Writer_Helper<1, In_Data_Type>

/// @details @em Specialization_Fix_Len_String: @p Mem_Type_Class = 2.
#ifndef DOXY
template <typename In_Data_Type>
struct Writer_Helper<2, In_Data_Type>
#else
template <int Mem_Type_Class, typename In_Data_Type>
struct Writer_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Fix_Len_String.
     * Write as fixed-length or variable-length strings.
     * @param grp_id Parent group object.
     * @param loc_name dataset/attribute name.
     * @param as_ds Flag; create dataset iff true.
     * @param dspace_id HDF5 dataspace.
     * @param sz Number of elements.
     * @param in Source memory address.
     * @param file_dtype_id HDF5 file datatype:
     * if 0, use @p mem_dtype_id;
     * if >0, write as fixed-length strings of this size;
     * if <0, write as variable-length strings.
     */
    void operator () (hid_t grp_id, std::string const & loc_name,
                      bool as_ds, hid_t dspace_id, size_t sz,
                      In_Data_Type const * in, hid_t file_dtype_id = 0) const
    {
        HDF_Object_Holder mem_dtype_id_holder;
        HDF_Object_Holder file_dtype_id_holder;
        std::vector<char const *> charptr_buff;
        void const * vptr_in = in;
        if (file_dtype_id >= 0)
        {
            mem_dtype_id_holder = Util::make_str_type(sizeof(In_Data_Type));
            if (file_dtype_id == 0)
            {
                file_dtype_id = mem_dtype_id_holder.id;
            }
            else // file_dtype_id > 0
            {
                file_dtype_id_holder = Util::make_str_type(file_dtype_id);
                file_dtype_id = file_dtype_id_holder.id;
            }
        }
        else // file_dtype_id < 0: write as varlen strings
        {
            mem_dtype_id_holder = Util::make_str_type(-1);
            file_dtype_id = mem_dtype_id_holder.id;
            // prepare array of pointers
            charptr_buff.resize(sz);
            for (hsize_t i = 0; i < sz; ++i)
            {
                charptr_buff[i] = &in[i][0];
            }
            vptr_in = charptr_buff.data();
        }
        Writer_Base::create_and_write(
            grp_id, loc_name, as_ds,
            dspace_id, mem_dtype_id_holder.id, file_dtype_id,
            vptr_in);
    }
}; // struct Writer_Helper<2, In_Data_Type>

/// @details @em Specialization_Var_Len_String: @p Mem_Type_Class = 3.
#ifndef DOXY
template <>
struct Writer_Helper<3, std::string>
#else
template <int Mem_Type_Class, typename In_Data_Type>
struct Writer_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Var_Len_String.
     * Write as fixed-length or variable-length strings.
     * @param grp_id Parent group object.
     * @param loc_name dataset/attribute name.
     * @param as_ds Flag; create dataset iff true.
     * @param dspace_id HDF5 dataspace.
     * @param sz Number of elements.
     * @param in Source memory address.
     * @param file_dtype_id HDF5 file datatype:
     * if -1, write variable-length strings;
     * if >=0, write as fixed-length strings of that size.
     */
    void operator () (hid_t grp_id, std::string const & loc_name,
                      bool as_ds, hid_t dspace_id, size_t sz,
                      std::string const * in, hid_t file_dtype_id = -1) const
    {
        HDF_Object_Holder mem_dtype_id_holder;
        std::vector<char const *> charptr_buff;
        std::vector<char> char_buff;
        void const * vptr_in;
        if (file_dtype_id == -1) // varlen to varlen
        {
            mem_dtype_id_holder = Util::make_str_type(-1);
            charptr_buff.resize(sz);
            for (hsize_t i = 0; i < sz; ++i)
            {
                charptr_buff[i] = in[i].data();
            }
            vptr_in = charptr_buff.data();
        }
        else // varlen to fixlen
        {
            assert(file_dtype_id > 0 or sz == 1); // file_dtype_id == 0 only allowed for single strings
            size_t slen = file_dtype_id > 0 ? file_dtype_id : in[0].size() + 1;
            assert(slen <= std::numeric_limits<long>::max());
            mem_dtype_id_holder = Util::make_str_type(slen);
            char_buff.resize(sz * slen);
            for (hsize_t i = 0; i < sz; ++i)
            {
                for (size_t j = 0; j < slen - 1; ++j)
                {
                    char_buff[i * slen + j] = j < in[i].size()? in[i][j] : '\0';
                }
                char_buff[i * slen + slen - 1] = '\0';
            }
            vptr_in = char_buff.data();
        }
        Writer_Base::create_and_write(
            grp_id, loc_name, as_ds,
            dspace_id, mem_dtype_id_holder.id, mem_dtype_id_holder.id,
            vptr_in);
    }
}; // struct Writer_Helper<3, std::string>

/// @details @em Specialization_Compound: @p Mem_Type_Class = 4.
#ifndef DOXY
template <typename In_Data_Type>
struct Writer_Helper<4, In_Data_Type>
#else
template <int Mem_Type_Class, typename In_Data_Type>
struct Writer_Helper
#endif
{
    /**
     * Functor operator. @em Specialization_Compound.
     * Go through the list of compound members;
     * write all of the ones that don't need conversion in one go;
     * write the remaining ones one at a time.
     * @param grp_id Parent group object.
     * @param loc_name dataset/attribute name.
     * @param as_ds Flag; create dataset iff true.
     * @param dspace_id HDF5 dataspace.
     * @param sz Number of elements.
     * @param in Source memory address.
     * @param cm Compound map.
     */
    void operator () (hid_t grp_id, std::string const & loc_name,
                      bool as_ds, hid_t dspace_id, size_t sz,
                      In_Data_Type const * in, Compound_Map const & cm) const
    {
        HDF_Object_Holder obj_id_holder;
        // create object
        {
            // create the file type
            HDF_Object_Holder file_dtype_id_holder(
                cm.build_type(sizeof(In_Data_Type), nullptr, false));
            obj_id_holder = Writer_Base::create(
                grp_id, loc_name, as_ds,
                dspace_id, file_dtype_id_holder.id);
        }
        // define functor that selects members which can be written with implicit conversion
        auto implicit_conversion = [] (detail::Compound_Member_Description const & e) {
            return (e.is_numeric() or e.is_char_array());
        };
        // write fields which do not need conversion, all-in-one
        {
            HDF_Object_Holder mem_dtype_id_holder(
                cm.build_type(sizeof(In_Data_Type), implicit_conversion, true));
            Writer_Base::write(obj_id_holder.id, as_ds, mem_dtype_id_holder.id, in);
        }
        // write fields which need conversion, one-by-one
        {
            auto mptr_l = cm.get_member_ptr_list();
            for (auto const & p : mptr_l)
            {
                auto const & e = *p.first.back();
                if (implicit_conversion(e)) continue;
                if (not as_ds) throw Exception("string in compound is supported in datasets, but not attributes");
                size_t mem_offset = p.second;
                if (e.is_string())
                {
                    // prepare memory vector of char*
                    std::vector<char const *> charptr_buff(sz);
                    for (size_t i = 0; i < sz; ++i)
                    {
                        charptr_buff[i] = reinterpret_cast<std::string const *>(
                            reinterpret_cast<char const *>(&in[i]) + mem_offset)->data();
                    }
                    // create flat hdf5 type
                    //HDF_Object_Holder mem_dtype_id_holder(Compound_Map::build_flat_type(p.first));
                    HDF_Object_Holder mem_dtype_id_holder(
                        cm.build_type(sizeof(In_Data_Type),
                                      [&e] (detail::Compound_Member_Description const & _e) {
                                          return &_e == &e;
                                      },
                                      false));
                    Writer_Base::write(obj_id_holder.id, as_ds, mem_dtype_id_holder.id, charptr_buff.data());
                }
            }
        }
    }
}; // struct Writer_Helper<4, In_Data_Type>

/**
 * Prepare for a write operation.
 * Branch statically on the source data type @p Src_Type:
 * if object is a vector, write a simple extent;
 * if object is not a vector, write a scalar.
 * Pass control to appropriate @p Writer_Helper.
 * Note: called by @p File::write.
 */
template <typename Src_Type>
struct Writer;

/**
 * @struct Writer
 * @details @em Specialization_Default: @p Src_Type = @p In_Data_Type.
 */
#ifndef DOXY
template <typename In_Data_Type>
struct Writer
#else
template <typename Src_Type>
struct Writer
#endif
{
    /**
     * Functor operator.
     * @em Specialization_Default.
     * Create dataspace and invoke @p Writer_Helper.
     * @param grp_id HDF5 parent group
     * @param loc_name Object name to write
     * @param as_ds Flag; create dataset iff true.
     * @param in Source (single address).
     * @param args Optional reading arguments passed to @p Writer_Helper.
     */
    template <typename ...Args>
    void operator () (hid_t grp_id, std::string const & loc_name,
                      bool as_ds, In_Data_Type const & in, Args && ...args) const
    {
        // create dataspace
        HDF_Object_Holder dspace_id_holder(
            Util::wrap(H5Screate, H5S_SCALAR),
            Util::wrapped_closer(H5Sclose));
        Writer_Helper<mem_type_class<In_Data_Type>::value, In_Data_Type>()(
            grp_id, loc_name, as_ds,
            dspace_id_holder.id, 1,
            &in, std::forward<Args>(args)...);
    }
}; // struct Writer<In_Data_Type>
/**
 * @struct Writer
 * @details @em Specialization_Vector: @p Src_Type = @p std::vector<In_Data_Type>.
 */
#ifndef DOXY
template <typename In_Data_Type>
struct Writer<std::vector<In_Data_Type>>
#else
template <typename Src_Type>
struct Writer
#endif
{
    /**
     * Functor operator.
     * @em Specialization_Vector.
     * Create dataspace and invoke @p Writer_Helper.
     * @param grp_id HDF5 parent group
     * @param loc_name Object name to write
     * @param as_ds Flag; create dataset iff true.
     * @param in Source.
     * @param args Optional reading arguments passed to @p Writer_Helper.
     */
    template <typename ...Args>
    void operator () (hid_t grp_id, std::string const & loc_name,
                      bool as_ds, std::vector<In_Data_Type> const & in, Args && ...args) const
    {
        assert(not in.empty());
        // create dataspace
        hsize_t sz = in.size();
        HDF_Object_Holder dspace_id_holder(
            Util::wrap(H5Screate_simple, 1, &sz, nullptr),
            Util::wrapped_closer(H5Sclose));
        Writer_Helper<mem_type_class<In_Data_Type>::value, In_Data_Type>()(
            grp_id, loc_name, as_ds,
            dspace_id_holder.id, sz,
            in.data(), std::forward<Args>(args)...);
    }
}; // struct Writer<std::vector<In_Data_Type>>

} // namespace detail

/// An HDF5 file reader
class File
{
public:
    /// A map of attributes of string type.
    typedef std::map<std::string, std::string> Attr_Map;

    /// Ctor: default
    File() : _file_id(0) {}
    /**
     * Ctor: from file name
     * @param file_name File name to open.
     * @param rw Flag: open for writing iff true.
     */
    File(std::string const & file_name, bool rw = false) : _file_id(0) { open(file_name, rw); }
    /// Ctor: copy
    File(File const &) = delete;
    /// Asop: copy
    File & operator = (File const &) = delete;
    /// Dtor
    ~File() { if (is_open()) close(); }

    /// Check if file is open.
    bool is_open() const { return _file_id > 0; }
    /// Check if file is open for writing.
    bool is_rw() const { return _rw; }
    /// Get file name.
    std::string const & file_name() const { return _file_name; }

    /**
     * Create file.
     * @param file_name File name to create.
     * @param truncate Control behaviour if file exists: if true, truncate; if false, fail.
     */
    void create(std::string const & file_name, bool truncate = false)
    {
        if (is_open()) close();
        _file_name = file_name;
        _rw = true;
        _file_id = H5Fcreate(file_name.c_str(), truncate? H5F_ACC_TRUNC : H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
        if (not is_open()) throw Exception(_file_name + ": error in H5Fcreate");
    } // create()
    /**
     * Open file.
     * @param file_name File name to open.
     * @param rw Flag: open for writing iff true.
     */
    void open(std::string const & file_name, bool rw = false)
    {
        if (is_open()) close();
        _file_name = file_name;
        _rw = rw;
        _file_id = H5Fopen(file_name.c_str(), not rw? H5F_ACC_RDONLY : H5F_ACC_RDWR, H5P_DEFAULT);
        if (not is_open()) throw Exception(_file_name + ": error in H5Fopen");
    } // open()
    /// Close file
    void close()
    {
        if (not is_open()) return;
        if (H5Fget_obj_count(_file_id, H5F_OBJ_ALL | H5F_OBJ_LOCAL) != 1) throw Exception(_file_name + ": HDF5 memory leak");
        int status = H5Fclose(_file_id);
        if (status < 0) throw Exception(_file_name + ": error in H5Fclose");
        _file_id = 0;
        _file_name.clear();
    } // close()
    /// Check if file name is a valid HDF5 file.
    static bool is_valid_file(std::string const & file_name)
    {
        std::ifstream ifs(file_name);
        if (not ifs) return false;
        (void)ifs.peek();
        if (not ifs) return false;
        ifs.close();
        auto status = H5Fis_hdf5(file_name.c_str());
        if (status <= 0) return 0;
        auto file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); // error if file is truncated
        if (file_id < 0) return 0;
        status = H5Fclose(file_id);
        if (status < 0) throw Exception(file_name + ": error in H5Fclose");
        return 1;
    } // is_valid_file()

    /// Get HDF5 object count.
    static int get_object_count()
    {
        return H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
    } // get_object_count()

    /**
     * Check if an object exists that is a group.
     * @param loc_full_name Full path.
     */
    bool
    group_exists(std::string const & loc_full_name) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name.front() == '/');
        if (loc_full_name == "/") return true;
        auto && loc = split_full_name(loc_full_name);
        // check all path elements exist, except for what is to the right of the last '/'
        // sets active path
        return path_exists(loc.first) and check_object_type(loc_full_name, H5O_TYPE_GROUP);
    } // group_exists()
    /**
     * Check if an object exists that is a dataset.
     * @param loc_full_name Full path.
     */
    bool
    dataset_exists(std::string const & loc_full_name) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name.front() == '/');
        if (loc_full_name == "/") return false;
        auto && loc = split_full_name(loc_full_name);
        // check all path elements exist, except for what is to the right of the last '/'
        // sets active path
        return path_exists(loc.first) and check_object_type(loc_full_name, H5O_TYPE_DATASET);
    } // dataset_exists()
    /**
     * Check if an object exists that is a group or a dataset.
     * @param loc_full_name Full path.
     */
    bool
    group_or_dataset_exists(std::string const & loc_full_name) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name.front() == '/');
        if (loc_full_name == "/") return true;
        auto && loc = split_full_name(loc_full_name);
        // check all path elements exist, except for what is to the right of the last '/'
        // sets active path
        return (path_exists(loc.first) and
                (check_object_type(loc_full_name, H5O_TYPE_DATASET) or
                 check_object_type(loc_full_name, H5O_TYPE_GROUP)));
    } // group_or_dataset_exists()
    /**
     * Check if an object exists that is an attribute.
     * @param loc_full_name Full path.
     */
    bool
    attribute_exists(std::string const & loc_full_name) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name.front() == '/');
        if (loc_full_name == "/") return false;
        auto && loc = split_full_name(loc_full_name);
        // check all path elements exist, except for what is to the right of the last '/'
        // sets active path
        if (not group_or_dataset_exists(loc.first)) return false;
        // check if target is an attribute
        int status = H5Aexists_by_name(_file_id, loc.first.c_str(), loc.second.c_str(), H5P_DEFAULT);
        if (status < 0) throw Exception("error in H5Aexists_by_name");
        return status > 0;
    } // attribute_exists()
    /**
     * Check if an object exists.
     * @param loc_full_name Full path.
     */
    bool
    exists(std::string const & loc_full_name) const
    {
        return attribute_exists(loc_full_name) or dataset_exists(loc_full_name);
    } // exists()

    /**
     * Read dataset/attribute.
     * @param loc_full_name Full path.
     * @param out Destination (single address or vector reference).
     * @param args Extra arguments to pass on to @p Reader.
     */
    template <typename Data_Storage, typename ...Args>
    void
    read(std::string const & loc_full_name, Data_Storage & out, Args && ...args) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name[0] == '/');
        auto && loc = split_full_name(loc_full_name);
        Exception::active_path() = loc_full_name;
        detail::HDF_Object_Holder grp_id_holder(
            detail::Util::wrap(H5Oopen, _file_id, loc.first.c_str(), H5P_DEFAULT),
            detail::Util::wrapped_closer(H5Oclose));
        detail::Reader<Data_Storage>()(grp_id_holder.id, loc.second,
                                         out, std::forward<Args>(args)...);
    } // read()
    /**
     * Write dataset/attribute.
     * @param loc_full_name Full path.
     * @param as_ds Flag; create dataset iff true.
     * @param in Source (single address or vector reference).
     * @param args Extra arguments to pass on to @p Writer.
     */
    template <typename In_Data_Storage, typename ...Args>
    void
    write(std::string const & loc_full_name, bool as_ds, In_Data_Storage const & in, Args && ...args) const
    {
        assert(is_open());
        assert(is_rw());
        assert(not loc_full_name.empty() and loc_full_name[0] == '/');
        assert(not exists(loc_full_name));
        auto && loc = split_full_name(loc_full_name);
        Exception::active_path() = loc_full_name;
        detail::HDF_Object_Holder grp_id_holder;
        if (group_or_dataset_exists(loc.first))
        {
            grp_id_holder = detail::HDF_Object_Holder(
                detail::Util::wrap(H5Oopen, _file_id, loc.first.c_str(), H5P_DEFAULT),
                detail::Util::wrapped_closer(H5Oclose));
        }
        else
        {
            detail::HDF_Object_Holder lcpl_id_holder(
                detail::Util::wrap(H5Pcreate, H5P_LINK_CREATE),
                detail::Util::wrapped_closer(H5Pclose));
            detail::Util::wrap(H5Pset_create_intermediate_group, lcpl_id_holder.id, 1);
            grp_id_holder = detail::HDF_Object_Holder(
                detail::Util::wrap(H5Gcreate2, _file_id, loc.first.c_str(), lcpl_id_holder.id, H5P_DEFAULT, H5P_DEFAULT),
                detail::Util::wrapped_closer(H5Gclose));
        }
        detail::Writer<In_Data_Storage>()(grp_id_holder.id, loc.second, as_ds, in, std::forward<Args>(args)...);
    } // write()
    /**
     * Write dataset.
     * @param loc_full_name Full path.
     * @param in Source (single address or vector reference).
     * @param args Extra arguments to pass on to @p Writer.
     */
    template <typename In_Data_Storage, typename ...Args>
    void
    write_dataset(std::string const & loc_full_name, In_Data_Storage const & in, Args && ...args) const
    {
        write(loc_full_name, true, in, std::forward<Args>(args)...);
    } // write_dataset()
    /**
     * Write attribute.
     * @param loc_full_name Full path.
     * @param in Source (single address or vector reference).
     * @param args Extra arguments to pass on to @p Writer.
     */
    template <typename In_Data_Storage, typename ...Args>
    void
    write_attribute(std::string const & loc_full_name, In_Data_Storage const & in, Args && ...args) const
    {
        write(loc_full_name, false, in, std::forward<Args>(args)...);
    } // write_attribute()

    /**
     * List group.
     * Return a list of names (groups/datasets) in the given group.
     * @param group_full_name Full path.
     */
    std::vector<std::string>
    list_group(std::string const & group_full_name) const
    {
        std::vector<std::string> res;
        Exception::active_path() = group_full_name;
        assert(group_exists(group_full_name));
        detail::HDF_Object_Holder g_id_holder(
            detail::Util::wrap(H5Gopen2, _file_id, group_full_name.c_str(), H5P_DEFAULT),
            detail::Util::wrapped_closer(H5Gclose));
        H5G_info_t g_info;
        detail::Util::wrap(H5Gget_info, g_id_holder.id, &g_info);
        res.resize(g_info.nlinks);
        for (unsigned i = 0; i < res.size(); ++i)
        {
            // find size first
            long sz1 = detail::Util::wrap(H5Lget_name_by_idx, _file_id, group_full_name.c_str(),
                                          H5_INDEX_NAME, H5_ITER_NATIVE, i, nullptr, 0, H5P_DEFAULT);
            res[i].resize(sz1);
            long sz2 = detail::Util::wrap(H5Lget_name_by_idx, _file_id, group_full_name.c_str(),
                                          H5_INDEX_NAME, H5_ITER_NATIVE, i, &res[i][0], sz1+1, H5P_DEFAULT);
            if (sz1 != sz2) throw Exception("error in H5Lget_name_by_idx: sz1!=sz2");
        }
        return res;
    } // list_group()
    /**
     * List attributes.
     * Return a list of attribute names of the given object.
     * @param loc_full_name Full path.
     */
    std::vector<std::string>
    get_attr_list(std::string const & loc_full_name) const
    {
        std::vector<std::string> res;
        Exception::active_path() = loc_full_name;
        assert(group_or_dataset_exists(loc_full_name));
        detail::HDF_Object_Holder id_holder(
            detail::Util::wrap(H5Oopen, _file_id, loc_full_name.c_str(), H5P_DEFAULT),
            detail::Util::wrapped_closer(H5Oclose));
        H5O_info_t info;
        detail::Util::wrap(H5Oget_info, id_holder.id, &info);
        // num_attrs in info.num_attrs
        for (unsigned i = 0; i < (unsigned)info.num_attrs; ++i)
        {
            int name_sz = detail::Util::wrap(H5Aget_name_by_idx, id_holder.id, ".",
                                             H5_INDEX_NAME, H5_ITER_NATIVE, i, nullptr, 0, H5P_DEFAULT);
            std::string tmp(name_sz, '\0');
            detail::Util::wrap(H5Aget_name_by_idx, id_holder.id, ".",
                               H5_INDEX_NAME, H5_ITER_NATIVE, i, &tmp[0], name_sz + 1, H5P_DEFAULT);
            res.emplace_back(std::move(tmp));
        }
        return res;
    } // get_attr_list()
    /**
     * Read attribute map.
     * @param path Full path.
     * @param recurse Flag: if true, recurse into subgroups.
     */
    Attr_Map
    get_attr_map(std::string const & path, bool recurse = false) const
    {
        Attr_Map res;
        std::queue<std::string> q;
        q.push("");
        while (not q.empty())
        {
            auto pt = q.front();
            q.pop();
            auto full_path = pt.empty()? path : path + "/" + pt;
            auto a_list = get_attr_list(full_path);
            for (auto const & a : a_list)
            {
                std::string tmp;
                read(full_path + "/" + a, tmp);
                res[pt.empty()? a : pt + "/" + a] = tmp;
            }
            if (recurse and group_exists(full_path))
            {
                auto sg_l = list_group(full_path);
                for (auto const & sg : sg_l)
                {
                    q.push(pt.empty()? sg : pt + "/" + sg);
                }
            }
        }
        return res;
    } // get_attr_map()
    /**
     * Write attribute map.
     * @param path Full path.
     * @param attr_m Attribute map.
     */
    void
    add_attr_map(std::string const & path, Attr_Map const & attr_m) const
    {
        for (auto const & p : attr_m)
        {
            write_attribute(path + "/" + p.first, p.second);
        }
    } // add_attr_map()
    /**
     * List field names.
     * Return a list of struct field names in the given compound dataset/attribute.
     * @param loc_full_name Full path.
     */
    std::vector<std::string>
    get_struct_members(std::string const & loc_full_name) const
    {
        std::vector<std::string> res;
        Exception::active_path() = loc_full_name;
        assert(attribute_exists(loc_full_name) or dataset_exists(loc_full_name));
        detail::HDF_Object_Holder attr_id_holder;
        detail::HDF_Object_Holder ds_id_holder;
        detail::HDF_Object_Holder type_id_holder;
        if (attribute_exists(loc_full_name))
        {
            auto && loc = split_full_name(loc_full_name);
            attr_id_holder = detail::HDF_Object_Holder(
                detail::Util::wrap(H5Aopen_by_name, _file_id, loc.first.c_str(), loc.second.c_str(),
                                   H5P_DEFAULT, H5P_DEFAULT),
                detail::Util::wrapped_closer(H5Aclose));
            type_id_holder = detail::HDF_Object_Holder(
                detail::Util::wrap(H5Aget_type, attr_id_holder.id),
                detail::Util::wrapped_closer(H5Tclose));
        }
        else
        {
            ds_id_holder = detail::HDF_Object_Holder(
                detail::Util::wrap(H5Oopen, _file_id, loc_full_name.c_str(), H5P_DEFAULT),
                detail::Util::wrapped_closer(H5Oclose));
            type_id_holder = detail::HDF_Object_Holder(
                detail::Util::wrap(H5Dget_type, ds_id_holder.id),
                detail::Util::wrapped_closer(H5Tclose));
        }
        if (detail::Util::wrap(H5Tget_class, type_id_holder.id) == H5T_COMPOUND)
        {
            // type is indeed a struct
            int nmem = detail::Util::wrap(H5Tget_nmembers, type_id_holder.id);
            for (int i = 0; i < nmem; ++i)
            {
                char* s = detail::Util::wrap(H5Tget_member_name, type_id_holder.id, i);
                res.emplace_back(s);
                free(s);
            }
        }
        return res;
    } // get_struct_members()

    /**
     * Copy an attribute between files.
     * @param src_f Source file.
     * @param dst_f Destination file.
     * @param src_full_path Source path.
     * @param _dst_full_path Destination path; if empty, use @p src_full_path.
     */
    static void
    copy_attribute(File const & src_f, File const & dst_f,
                   std::string const & src_full_path, std::string const & _dst_full_path = std::string())
    {
        if (not src_f.is_open()) throw Exception("source file not open");
        if (not dst_f.is_open()) throw Exception("destination file not open");
        if (not dst_f.is_rw()) throw Exception("destination file not writeable");
        std::string const & dst_full_path = (_dst_full_path.empty()? src_full_path : _dst_full_path);
        if (not src_f.attribute_exists(src_full_path)) throw Exception("source attribute missing");
        if (dst_f.group_or_dataset_exists(dst_full_path) or
            dst_f.attribute_exists(dst_full_path)) throw Exception("destination path exists");
        // compute paths
        auto && src_path = split_full_name(src_full_path);
        auto && dst_path = split_full_name(dst_full_path);
        // open source attribute
        detail::HDF_Object_Holder src_attr_id_holder(
            detail::Util::wrap(H5Aopen_by_name, src_f._file_id, src_path.first.c_str(), src_path.second.c_str(),
                               H5P_DEFAULT, H5P_DEFAULT),
            detail::Util::wrapped_closer(H5Aclose));
        // open source attribute datatype
        detail::HDF_Object_Holder src_attr_dtype_id_holder(
            detail::Util::wrap(H5Aget_type, src_attr_id_holder.id),
            detail::Util::wrapped_closer(H5Tclose));
        if (H5Tget_class(src_attr_dtype_id_holder.id) == H5T_INTEGER)
        {
            if (H5Tget_sign(src_attr_dtype_id_holder.id) == H5T_SGN_NONE)
            {
                unsigned long long tmp;
                src_f.read(src_full_path, tmp);
                dst_f.write_attribute(dst_full_path, tmp, src_attr_dtype_id_holder.id);
            }
            else if (H5Tget_sign(src_attr_dtype_id_holder.id) == H5T_SGN_2)
            {
                long long tmp;
                src_f.read(src_full_path, tmp);
                dst_f.write_attribute(dst_full_path, tmp, src_attr_dtype_id_holder.id);
            }
            else
            {
                throw Exception("error in H5Tget_sign");
            }
        }
        else if (H5Tget_class(src_attr_dtype_id_holder.id) == H5T_FLOAT)
        {
            long double tmp;
            src_f.read(src_full_path, tmp);
            dst_f.write_attribute(dst_full_path, tmp, src_attr_dtype_id_holder.id);
        }
        else if (H5Tget_class(src_attr_dtype_id_holder.id) == H5T_STRING)
        {
            std::string tmp;
            src_f.read(src_full_path, tmp);
            auto is_varlen = H5Tis_variable_str(src_attr_dtype_id_holder.id);
            if (is_varlen < 0) throw Exception("error in H5Tis_variable_str");
            if (is_varlen)
            {
                dst_f.write_attribute(dst_full_path, tmp, -1);
            }
            else
            {
                // not varlen; now deal with array-of-size-1 chars
                int sz = H5Tget_size(src_attr_dtype_id_holder.id);
                if (sz == 0) throw Exception("error in H5Tget_size");
                detail::HDF_Object_Holder src_attr_dspace_id_holder(
                    detail::Util::wrap(H5Aget_space, src_attr_id_holder.id),
                    detail::Util::wrapped_closer(H5Sclose));
                auto dspace_type = H5Sget_simple_extent_type(src_attr_dspace_id_holder.id);
                if (dspace_type == H5S_SCALAR)
                {
                    dst_f.write_attribute(dst_full_path, tmp, 0);
                }
                else if (dspace_type == H5S_SIMPLE)
                {
                    if (sz != 1) throw Exception("unsupported attribute type for copying: extent of string of size > 1");
                    std::vector<std::array<char, 1>> tmp_v(tmp.size());
                    for (unsigned i = 0; i < tmp.size(); ++i)
                    {
                        tmp_v[i][0] = tmp[i];
                    }
                    dst_f.write_attribute(dst_full_path, tmp_v);
                }
                else
                {
                    throw Exception("error in H5Sget_simple_extent_type");
                }
            }
        }
        else
        {
            throw Exception("unsupported attribute type for copying");
        }
    } // copy_attribute()

    /**
     * Copy attributes between files.
     * @param src_f Source file.
     * @param dst_f Destination file.
     * @param path Source path.
     * @param recurse Flag; if true, recurse.
     */
    static void
    copy_attributes(File const & src_f, File const & dst_f, std::string const & path, bool recurse = false)
    {
        auto a_l = src_f.get_attr_list(not path.empty()? path : std::string("/"));
        for (auto const & a : a_l)
        {
            copy_attribute(src_f, dst_f, path + "/" + a);
        }
        if (not recurse) return;
        auto sg_l = src_f.list_group(not path.empty()? path : std::string("/"));
        for (auto const & sg : sg_l)
        {
            if (src_f.group_exists(path + "/" + sg))
            {
                copy_attributes(src_f, dst_f, path + "/" + sg, true);
            }
        }
    } // copy_attributes()
private:
    std::string _file_name;
    hid_t _file_id;
    bool _rw;

    /**
     * Split a full name into path and name.
     * Note: @p full_name must begin with '/', and not end with '/' unless it equals "/".
     */
    static std::pair<std::string, std::string>
    split_full_name(std::string const & full_name)
    {
        assert(not full_name.empty() and
               full_name.front() == '/' and
               (full_name.size() == 1 or full_name.back() != '/'));
        if (full_name == "/") return std::make_pair(std::string("/"), std::string());
        auto pos = full_name.find_last_of('/');
        return (pos != std::string::npos
                ? std::make_pair(full_name.substr(0, pos > 0? pos : 1), full_name.substr(pos + 1))
                : std::make_pair(std::string(), std::string()));
    } // split_full_name()

    /// Determine if a path to an element exists
    bool
    path_exists(std::string const & full_path_name) const
    {
        assert(is_open());
        assert(not full_path_name.empty() and full_path_name.front() == '/');
        if (full_path_name == "/") return true;
        Exception::active_path() = full_path_name;
        // check all path elements exist, except for what is to the right of the last '/'
        size_t pos = 0;
        while (pos != std::string::npos)
        {
            ++pos;
            pos = full_path_name.find('/', pos);
            std::string tmp = full_path_name.substr(0, pos);
            // check link exists
            if (not detail::Util::wrap(H5Lexists, _file_id, tmp.c_str(), H5P_DEFAULT)) return false;
            // check object exists
            if (not detail::Util::wrap(H5Oexists_by_name, _file_id, tmp.c_str(), H5P_DEFAULT)) return false;
            // open object in order to check type
            detail::HDF_Object_Holder o_id_holder(
                detail::Util::wrap(H5Oopen, _file_id, tmp.c_str(), H5P_DEFAULT),
                detail::Util::wrapped_closer(H5Oclose));
            // check object is a group
            H5O_info_t o_info;
            detail::Util::wrap(H5Oget_info, o_id_holder.id, &o_info);
            if (o_info.type != H5O_TYPE_GROUP) return false;
        }
        return true;
    } // path_exists()

    /// Check if HDF5 object has given type
    bool
    check_object_type(std::string const & loc_full_name, H5O_type_t type_id) const
    {
        // check link exists
        if (loc_full_name != "/"
            and not detail::Util::wrap(H5Lexists, _file_id, loc_full_name.c_str(), H5P_DEFAULT)) return false;
        // check object exists
        if (not detail::Util::wrap(H5Oexists_by_name, _file_id, loc_full_name.c_str(), H5P_DEFAULT)) return false;
        // open object in order to check type
        detail::HDF_Object_Holder o_id_holder(
            detail::Util::wrap(H5Oopen, _file_id, loc_full_name.c_str(), H5P_DEFAULT),
            detail::Util::wrapped_closer(H5Oclose));
        // check object is a group
        H5O_info_t o_info;
        detail::Util::wrap(H5Oget_info, o_id_holder.id, &o_info);
        return o_info.type == type_id;
    } // check_object_type()
}; // class File

} // namespace hdf5_tools

#endif
