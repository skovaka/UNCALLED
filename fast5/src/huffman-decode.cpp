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
    string l;
    map< string, string > cw_v_id;
    vector< uint8_t > cw_v;
    while (getline(cin, l))
    {
        if (l[0] == '#')
        {
            istringstream iss(l.substr(1));
            string tmp0;
            string tmp1;
            getline(iss, tmp0, '=');
            iss >> tmp1;
            cw_v_id[tmp0] = tmp1;
        }
        else
        {
            unsigned x;
            istringstream(l) >> x;
            cw_v.push_back(x);
        }
    }
    auto val_v = hc.decode<int16_t>(cw_v, cw_v_id);
    for (auto x : val_v)
    {
        cout << x << endl;
    }
}
