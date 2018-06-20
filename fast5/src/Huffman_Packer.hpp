//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#ifndef __HUFFMAN_PACKER_HPP
#define __HUFFMAN_PACKER_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <stdexcept>
#include <cassert>
#include <bitset>

#include "logger.hpp"

namespace fast5
{

class Huffman_Packer
{
public:
    typedef std::vector< std::uint8_t > Code_Type;
    typedef std::map< std::string, std::string > Code_Params_Type;

    Huffman_Packer() = default;
    Huffman_Packer(Huffman_Packer const &) = delete;
    Huffman_Packer(Huffman_Packer &&) = default;
    Huffman_Packer & operator = (Huffman_Packer const &) = delete;
    Huffman_Packer & operator = (Huffman_Packer &&) = default;
    Huffman_Packer(std::istream & is, std::string const & cwm_name)
    {
        load_codeword_map(is, cwm_name);
    }
    Huffman_Packer(std::vector< std::string > const & v, std::string const & cwm_name)
    {
        load_codeword_map(v.begin(), v.end(), cwm_name);
    }
    template < typename Iterator >
    Huffman_Packer(Iterator it_begin, Iterator it_end, std::string const & cwm_name)
    {
        load_codeword_map(it_begin, it_end, cwm_name);
    }

    void load_codeword_map(std::istream & is, std::string const & cwm_name)
    {
        _cwm_name = cwm_name;
        std::string v_s;
        std::string cw_s;
        while (is >> v_s >> cw_s)
        {
            add_codeword(v_s, cw_s);
        }
    }
    template < typename Iterator >
    void load_codeword_map(Iterator it_begin, Iterator it_end, std::string const & cwm_name)
    {
        _cwm_name = cwm_name;
        for (auto it = it_begin; it != it_end and std::next(it) != it_end; it += 2)
        {
            add_codeword(*it, *next(it));
        }
    }

    template < typename Int_Type >
    std::pair< Code_Type, Code_Params_Type >
    encode(std::vector< Int_Type > const & v, bool encode_diff = false) const
    {
        Code_Type res;
        Code_Params_Type res_params = id();
        res_params["code_diff"] = encode_diff? "1" : "0";
        std::ostringstream oss;
        oss << v.size();
        res_params["size"] = oss.str();
        uint64_t buff = 0;
        uint8_t buff_len = 0;
        bool reset = true;
        Int_Type last = 0;
        unsigned i = 0;
        long long int val;
        long long int x;
        while (true)
        {
            assert(buff_len <= 64);
            // flush buffer
            while (buff_len >= 8)
            {
                res.push_back(buff & 0xFF);
                buff >>= 8;
                buff_len -= 8;
            }
            assert(buff_len < 8);
            if (reset)
            {
                assert(buff_len == 0);
                if (i == v.size()) break;
                //LOG(debug) << "absolute value val=" << v[i] << std::endl;
                for (unsigned j = 0; j < sizeof(Int_Type); ++j)
                {
                    std::uint8_t y = (v[i] >> (8 * j)) & 0xFF;
                    //LOG(debug) << "byte " << j << ": " << std::bitset<8>(y) << std::endl;
                    res.push_back(y);
                }
                reset = false;
                last = v[i];
                ++i;
            }
            else // not reset
            {
                if (i < v.size())
                {
                    val = v[i];
                    x = encode_diff? val - last : val;
                    reset = _cwm.count(x) == 0;
                    //LOG(debug) << "relative value: val=" << v[i] << " last=" << last << " x=" << x << " reset=" << reset << std::endl;
                }
                else
                {
                    reset = true;
                    //LOG(debug) << "end: reset=1" << std::endl;
                }
                auto p = (not reset? _cwm.at(x) : _cwm.at(break_cw()));
                buff |= (p.first << buff_len);
                buff_len += p.second;
                if (not reset)
                {
                    last = v[i];
                    ++i;
                }
                else if ((buff_len % 8) > 0) // and reset
                {
                    buff_len += 8 - (buff_len % 8);
                }
                        
            }
        }
        oss.str("");
        oss << std::fixed << std::setprecision(2) << (double)(res.size() * 8) / v.size();
        res_params["avg_bits"] = oss.str();
        return std::make_pair(std::move(res), std::move(res_params));
    }

    template < typename Int_Type >
    std::vector< Int_Type >
    decode(Code_Type const & v, Code_Params_Type const & v_params) const
    {
        check_params(v_params);
        bool decode_diff = v_params.at("code_diff") == "1";
        std::vector< Int_Type > res;
        std::uint64_t buff = 0;
        std::uint8_t buff_len = 0;
        bool reset = true;
        Int_Type last = 0;
        unsigned i = 0;
        while (i < v.size() or buff_len > 0)
        {
            assert(buff_len <= 64);
            // fill buffer
            while (i < v.size() and buff_len <= 56)
            {
                uint64_t y = v[i];
                buff |= (y << buff_len);
                buff_len += 8;
                ++i;
            }
            assert(buff_len <= 64);
            if (reset)
            {
                assert((buff_len % 8) == 0);
                assert(buff_len / 8 >= sizeof(Int_Type));
                //LOG(debug) << "absolute value" << std::endl;
                Int_Type x = 0;
                for (unsigned j = 0; j < sizeof(Int_Type); ++j)
                {
                    std::uint64_t y = (buff & 0xFF);
                    //LOG(debug) << "byte " << j << ": " << std::bitset<8>(y) << std::endl;
                    x |= (y << (8 * j));
                    buff >>= 8;
                    buff_len -= 8;
                }
                //LOG(debug) << "got: val=" << x << std::endl;
                res.push_back(x);
                last = x;
                reset = false;
            }
            else // not reset
            {
                //LOG(debug) << "reading relative value" << std::endl;
                // TODO: faster decoding
                // currently, try all codewords one by one
                auto it = _cwm.begin();
                while (it != _cwm.end())
                {
                    if ((buff & ((1llu << it->second.second) - 1)) == it->second.first)
                    {
                        break;
                    }
                    ++it;
                }
                if (it == _cwm.end())
                {
                    LOG_THROW
                        << "codeword not found: buff=" << std::bitset<64>(buff);
                }
                auto x = it->first;
                auto p = it->second;
                assert(buff_len >= p.second);
                buff >>= p.second;
                buff_len -= p.second;
                if (x != break_cw())
                {
                    //LOG(debug) << "got: x=" << x << " last=" << last << " val=" << x + last << " cw_len=" << (int)p.second << std::endl;
                    if (decode_diff) x += last;
                    if (sizeof(Int_Type) < 8
                        and (x < (long long)std::numeric_limits< Int_Type >::min()
                             or x > (long long)std::numeric_limits< Int_Type >::max()))
                    {
                        LOG_THROW
                            << "overflow";
                    }
                    res.push_back(x);
                    last = x;
                }
                else
                {
                    //LOG(debug) << "got: break cw_len=" << (int)p.second << std::endl;
                    reset = true;
                    if ((buff_len % 8) > 0)
                    {
                        buff >>= (buff_len % 8);
                        buff_len -= (buff_len % 8);
                    }
                }
            }
        }
        return res;
    }

    //
    // static coder access
    //
    static Huffman_Packer const &
    get_coder(std::string const & cwm_name)
    {
        static_init();
        if (cwm_m().count(cwm_name) == 0)
        {
            LOG_THROW
                << "missing codeword map: " + cwm_name;
        }
        return cwm_m().at(cwm_name);
    }

private:
    std::map< long long int, std::pair< std::uint64_t, std::uint8_t > > _cwm;
    std::string _cwm_name;
    static long long int break_cw()
    {
        static long long int const _break_cw = std::numeric_limits< long long int >::min();
        return _break_cw;
    }
    Code_Params_Type id() const
    {
        Code_Params_Type res;
        res["packer"] = "huffman_packer";
        res["format_version"] = "2";
        res["codeword_map_name"] = _cwm_name;
        return res;
    }
    void check_params(Code_Params_Type const & params) const
    {
        auto _id = id();
        if (params.at("packer") != _id.at("packer")
            or params.at("format_version") != _id.at("format_version")
            or params.at("codeword_map_name") != _id.at("codeword_map_name"))
        {
            LOG_THROW
                << "decode id mismatch";
        }
    }
    void add_codeword(std::string const & v_s, std::string const & cw_s)
    {
        long long int v;
        if (v_s != ".")
        {
            std::istringstream(v_s) >> v;
        }
        else
        {
            v = break_cw();
        }
        std::uint64_t cw = 0;
        if (cw_s.size() > 57)
        {
            LOG_THROW
                << "codeword too long: " + v_s + " " + cw_s;
        }
        std::uint8_t cw_l = cw_s.size();
        for (int i = cw_s.size() - 1; i >= 0; --i)
        {
            cw <<= 1;
            cw |= (cw_s[i] == '1');
        }
        _cwm[v] = std::make_pair(cw, cw_l);
    }

    static std::map< std::string, Huffman_Packer > & cwm_m()
    {
        static std::map< std::string, Huffman_Packer > _cwm_m;
        return _cwm_m;
    }
    static void static_init()
    {
        static bool inited = false;
        if (inited) return;
        std::deque< std::deque< std::string > > dd;
        dd.push_back(
#include "cwmap.fast5_rw_1.inl"
            );
        dd.push_back(
#include "cwmap.fast5_ed_skip_1.inl"
            );
        dd.push_back(
#include "cwmap.fast5_ed_len_1.inl"
            );
        dd.push_back(
#include "cwmap.fast5_fq_bp_1.inl"
            );
        dd.push_back(
#include "cwmap.fast5_fq_qv_1.inl"
            );
        dd.push_back(
#include "cwmap.fast5_ev_rel_skip_1.inl"
            );
        dd.push_back(
#include "cwmap.fast5_ev_move_1.inl"
            );
        cwm_m().clear();
        for (auto & d : dd)
        {
            auto cwm_name = d.front();
            Huffman_Packer hc(d.begin() + 1, d.end(), cwm_name);
            cwm_m()[cwm_name] = std::move(hc);
        }
        inited = true;
    } // static_init()
}; // class Huffman_Packer

} // namespace fast5

#endif
