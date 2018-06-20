//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#ifndef __BIT_PACKER_HPP
#define __BIT_PACKER_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <stdexcept>
#include <cassert>

#include "logger.hpp"

namespace fast5
{

class Bit_Packer
{
public:
    typedef std::vector< std::uint8_t > Code_Type;
    typedef std::map< std::string, std::string > Code_Params_Type;

    template < typename Int_Type >
    std::pair< Code_Type, Code_Params_Type >
    encode(std::vector< Int_Type > const & v, unsigned num_bits) const
    {
        Code_Type res;
        Code_Params_Type res_params;
        res_params["packer"] = "bit_packer";
        num_bits = std::min(num_bits, (unsigned)sizeof(Int_Type) * 8);
        std::ostringstream oss;
        oss << num_bits;
        res_params["num_bits"] = oss.str();
        oss.str("");
        oss << v.size();
        res_params["size"] = oss.str();
        long long unsigned buff = 0;
        unsigned buff_len = 0;
        auto val_mask = (1llu << num_bits) - 1;
        for (unsigned i = 0; i < v.size(); ++i)
        {
            // flush out buff
            while (buff_len >= 8)
            {
                res.push_back(buff & 0xFF);
                buff >>= 8;
                buff_len -= 8;
            }
            assert(buff_len < 8);
            long long unsigned x = v[i];
            if (buff_len + num_bits <= 64)
            {
                buff |= (x & val_mask) << buff_len;
                buff_len += num_bits;
            }
            else
            {
                assert(num_bits > 56);
                buff |= (x & 0xFF) << buff_len;
                res.push_back(buff & 0xFF);
                buff >>= 8;
                x >>= 8;
                buff |= (x & (val_mask >> 8)) << buff_len;
                buff_len += num_bits - 8;
            }
        }
        while (buff_len >= 8)
        {
            res.push_back(buff & 0xFF);
            buff >>= 8;
            buff_len -= 8;
        }
        if (buff_len > 0)
        {
            res.push_back(buff & 0xFF);
        }
        return std::make_pair(std::move(res), std::move(res_params));
    } // encode()

    template < typename Int_Type >
    std::vector< Int_Type >
    decode(Code_Type const & v, Code_Params_Type const & v_params) const
    {
        std::vector< Int_Type > res;
        unsigned num_bits;
        size_t sz;
        std::istringstream(v_params.at("num_bits")) >> num_bits;
        std::istringstream(v_params.at("size")) >> sz;
        if (v.size() != (sz * num_bits) / 8 + ((sz * num_bits) % 8 > 0? 1 : 0))
        {
            LOG_THROW
                << "incorrect size: v_size=" << v.size();
        }
        long long unsigned buff = 0;
        unsigned buff_len = 0;
        unsigned j = 0;
        auto val_mask = (1llu << num_bits) - 1;
        for (unsigned i = 0; i < sz; ++i)
        {
            while (j < v.size() and buff_len <= 64 - 8)
            {
                buff |= ((long long unsigned)v.at(j) << buff_len);
                ++j;
                buff_len += 8;
            }
            long long unsigned x;
            if (buff_len >= num_bits)
            {
                x = buff & val_mask;
                buff >>= num_bits;
                buff_len -= num_bits;
            }
            else
            {
                // 56 < buff_len < num_bits
                x = buff & 0xFF;
                buff >>= 8;
                buff_len -= 8;
                buff |= (v.at(j) << buff_len);
                ++j;
                buff_len += 8;
                x |= ((buff & (val_mask >> 8)) << 8);
                buff >>= (num_bits - 8);
                buff_len -= num_bits - 8;
            }
            res.push_back(x);
        }
        return res;
    } // decode()

    //
    // static packer access
    //
    static Bit_Packer const &
    get_packer()
    {
        static Bit_Packer _packer;
        return _packer;
    }
}; // class Bit_Packer

} // namespace fast5

#endif
