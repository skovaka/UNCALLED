//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "fast5_pack.hpp"
#include "logger.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    logger::Logger::set_default_level(logger::level::debug);
    if (argc != 2)
    {
        cerr << "use: " << argv[0] << " <codeword_file>" << endl;
        exit(EXIT_FAILURE);
    }
    string cw_fn = argv[1];
    ifstream cw_f(cw_fn);
    fast5_pack::Huffman_Diff_Coder hc(cw_f, cw_fn);
    int16_t x;
    std::vector< int16_t > val_v;
    while (cin >> x)
    {
        val_v.push_back(x);
    }
    auto p = hc.encode(val_v);
    for (auto const & p2 : p.second)
    {
        cout << "#" << p2.first << "=" << p2.second << endl;
    }
    for (auto y : p.first)
    {
        cout << (int)y << endl;
    }
}
