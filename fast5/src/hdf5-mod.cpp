//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#include <cassert>
#include <iostream>
#include <string>

#include "hdf5_tools.hpp"

using namespace std;
using namespace hdf5;

struct B
{
    int val_1;
    array< char, 6 > val_2;
    string val_3;
    friend ostream & operator << (ostream & os, const B & b)
    {
        os << "(val_1=" << b.val_1
           << ",val_2=\"" << string(b.val_2.begin(), b.val_2.end())
           << "\",val_3=\"" << b.val_3 << "\")";
        return os;
    }
};

struct A
{
    int val_1;
    int val_1a;
    float val_2;
    char val_3[6];
    array< char, 6 > val_4;
    string val_5;
    B val_6;
    friend ostream & operator << (ostream & os, const A & a)
    {
        os << "(val_1=" << a.val_1
           << ",val_1a=" << a.val_1a
           << ",val_2=" << a.val_2
           << ",val_3=\"" << a.val_3
           << "\",val_4=\"" << string(a.val_4.begin(), a.val_4.end())
           << "\",val_5=\"" << a.val_5
           << "\",val_6=" << a.val_6 << ")";
        return os;
    }
};

struct B_string
{
    string val_1;
    string val_2;
    string val_3;
    friend ostream & operator << (ostream & os, const B_string & b)
    {
        os << "(val_1=\"" << b.val_1
           << "\",val_2=\"" << b.val_2
           << "\",val_3=\"" << b.val_3 << "\")";
        return os;
    }
};
struct A_string
{
    string val_1;
    string val_1a;
    string val_2;
    string val_3;
    string val_4;
    string val_5;
    B_string val_6;
    friend ostream & operator << (ostream & os, const A_string & a)
    {
        os << "(val_1=\"" << a.val_1
           << "\",val_1a=\"" << a.val_1a
           << "\",val_2=\"" << a.val_2
           << "\",val_3=\"" << a.val_3
           << "\",val_4=\"" << a.val_4
           << "\",val_5=\"" << a.val_5
           << "\",val_6=" << a.val_6 << ")";
        return os;
    }
};

struct B_char_array
{
    array< char, 6 > val_2;
    array< char, 6 > val_3;
    friend ostream & operator << (ostream & os, const B_char_array & b)
    {
        os << "(val_2=\"" << string(b.val_2.begin(), b.val_2.end())
           << "\",val_3=\"" << string(b.val_3.begin(), b.val_3.end()) << ")";
        return os;
    }
};

struct A_char_array
{
    array< char, 6 > val_3;
    array< char, 6 > val_4;
    array< char, 6 > val_5;
    B_char_array val_6;
    friend ostream & operator << (ostream & os, const A_char_array & a)
    {
        os
           << "(val_3=\"" << string(a.val_3.begin(), a.val_3.end())
           << "\",val_4=\"" << string(a.val_4.begin(), a.val_4.end())
           << "\",val_5=\"" << string(a.val_5.begin(), a.val_5.end())
           << "\",val_6=" << a.val_6 << ")";
        return os;
    }
};

int main(int argc, char* argv[])
{
    if (argc != 2 and argc != 3)
    {
        cerr << "use: " << argv[0] << " [-f] <fast5_file>" << endl;
        return EXIT_FAILURE;
    }
    bool force = string(argv[1]) == "-f";
    string file_name(argv[force? 2 : 1]);
    {
        hdf5_tools::File f;
        //
        // All fast5 operations are performed inside a try-catch block. This should
        // resist various hdf5 errors without leaking memory.
        //
        try
        {
            //
            // create file; without -f, fail if it exist
            //
            f.create(file_name, force);
            assert(f.is_open());
            assert(f.is_rw());
            //
            // write a /file_version to allow f5ls to work
            //
            string file_version("42");
            f.write("/file_version", false, file_version);
            int val_1 = 42;
            float val_2 = 3.14;
            char val_3[6] = "ACGTA";
            array< char, 6 > val_4 = { "AACCG" };
            string val_5("CCCGG");
            static_assert(hdf5_tools::detail::mem_type_class< void >::value == 0, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_1) >::value == 1, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_2) >::value == 1, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_3) >::value == 2, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_4) >::value == 2, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_5) >::value == 3, "");
            static_assert(hdf5_tools::detail::mem_type_class< std::true_type >::value == 4, "");
            //
            // write integer
            //
            f.write("/val_1", false, val_1);
            f.write("/val_1_as_64", false, val_1, H5T_STD_I64LE);
            f.write("/val_1_v", false, vector< int >(3, val_1));
            //
            // write float
            //
            f.write("/val_2", false, val_2);
            f.write("/val_2_as_64", false, val_2, H5T_IEEE_F64LE);
            f.write("/val_2_v", false, vector< float >(3, val_2));
            //
            // write fixlen string: char[]
            //
            f.write("/val_3", false, val_3);
            f.write("/val_3_as_len_3", false, val_3, 3);
            f.write("/val_3_as_varlen", false, val_3, -1);
            //
            // write fixlen string: std::array< char >
            //
            f.write("/val_4", false, val_4);
            f.write("/val_4_as_len_3", false, val_4, 3);
            f.write("/val_4_as_varlen", false, val_4, -1);
            f.write("/val_4_v", false, vector < decltype(val_4) >(3, val_4));
            f.write("/val_4_v_as_len_3", false, vector < decltype(val_4) >(3, val_4), 3);
            f.write("/val_4_v_as_varlen", false, vector < decltype(val_4) >(3, val_4), -1);
            //
            // write varlen string
            //
            f.write("/val_5", false, val_5);
            f.write("/val_5_as_len_3", false, val_5, 3);
            f.write("/val_5_as_fixlen", false, val_5, 0);
            f.write("/val_5_v", false, vector< decltype(val_5) >(3, val_5));
            f.write("/val_5_v_as_len_3", false, vector< decltype(val_5) >(3, val_5), 3);
            f.write("/val_5_v_as_fixlen", false, vector< decltype(val_5) >(1, val_5), 0); // only size 1
            //
            // write compound
            //
            A val_6{ 1, 2, 3.14, "ACGTA", "CGTAC", "CCGGT", { 42, "GTTAC", "TTATT" } };
            hdf5_tools::Compound_Map cm_A;
            hdf5_tools::Compound_Map cm_B;
            cm_B.add_member("val_1", &B::val_1);
            cm_B.add_member("val_2", &B::val_2);
            cm_B.add_member("val_3", &B::val_3);
            for (const auto& e : cm_B.members())
            {
                clog << "cm_B: (" << (void*)&e << ")" << e << endl;
            }
            cm_A.add_member("val_1", &A::val_1);
            cm_A.add_member("val_2", &A::val_2);
            cm_A.add_member("val_3", &A::val_3);
            cm_A.add_member("val_4", &A::val_4);
            cm_A.add_member("val_5", &A::val_5);
            cm_A.add_member("val_6", &A::val_6, &cm_B);
            for (const auto& e : cm_A.members())
            {
                clog << "cm_A: (" << (void*)&e << ")" << e << endl;
            }
            auto l = cm_A.get_member_ptr_list();
            for (const auto& p : l)
            {
                clog << "member:";
                for (const auto& e_ptr : p.first)
                {
                    clog << " " << *e_ptr;
                }
                clog << "; total_offset=" << p.second << endl;
            }
            //f.write("/val_6a", false, val_6, cm_A);
            f.write("/val_6d", true, val_6, cm_A);
            vector< A > src(3, val_6);
            f.write("/val_6d_v", true, src, cm_A);
            clog << "wrote val_6d_v:" << endl;
            for (const auto& a : src)
            {
                clog << a << endl;
            }

            // reopen for reading
            f.close();
            f.open(file_name);

            //
            // test reading compound
            //
            // using original map
            {
                std::vector< A > dest;
                f.read("/val_6d_v", dest, cm_A);
                clog << "read val_6d_v:" << endl;
                for (const auto& a : dest)
                {
                    clog << a << endl;
                }
            }
            // using all strings
            {
                hdf5_tools::Compound_Map cm_A_string;
                hdf5_tools::Compound_Map cm_B_string;
                cm_B_string.add_member("val_1", &B_string::val_1);
                cm_B_string.add_member("val_2", &B_string::val_2);
                cm_B_string.add_member("val_3", &B_string::val_3);
                for (const auto& e : cm_B_string.members())
                {
                    clog << "cm_B_string: (" << (void*)&e << ")" << e << endl;
                }
                cm_A_string.add_member("val_1", &A_string::val_1);
                cm_A_string.add_member("val_2", &A_string::val_2);
                cm_A_string.add_member("val_3", &A_string::val_3);
                cm_A_string.add_member("val_4", &A_string::val_4);
                cm_A_string.add_member("val_5", &A_string::val_5);
                cm_A_string.add_member("val_6", &A_string::val_6, &cm_B_string);
                for (const auto& e : cm_A_string.members())
                {
                    clog << "cm_A_string: (" << (void*)&e << ")" << e << endl;
                }
                std::vector< A_string > dest;
                f.read("/val_6d_v", dest, cm_A_string);
                clog << "read val_6d_v using all-strings:" << endl;
                for (const auto& a : dest)
                {
                    clog << a << endl;
                }
            }
            // using char arrays
            {
                hdf5_tools::Compound_Map cm_A_char_array;
                hdf5_tools::Compound_Map cm_B_char_array;
                cm_B_char_array.add_member("val_2", &B_char_array::val_2);
                cm_B_char_array.add_member("val_3", &B_char_array::val_3);
                for (const auto& e : cm_B_char_array.members())
                {
                    clog << "cm_B_char_array: (" << (void*)&e << ")" << e << endl;
                }
                cm_A_char_array.add_member("val_3", &A_char_array::val_3);
                cm_A_char_array.add_member("val_4", &A_char_array::val_4);
                cm_A_char_array.add_member("val_5", &A_char_array::val_5);
                cm_A_char_array.add_member("val_6", &A_char_array::val_6, &cm_B_char_array);
                for (const auto& e : cm_A_char_array.members())
                {
                    clog << "cm_A_char_array: (" << (void*)&e << ")" << e << endl;
                }
                std::vector< A_char_array > dest;
                f.read("/val_6d_v", dest, cm_A_char_array);
                clog << "read val_6d_v using char arrays:" << endl;
                for (const auto& a : dest)
                {
                    clog << a << endl;
                }

            }
            f.close();
        }
        catch (hdf5_tools::Exception& e)
        {
            cout << "hdf5 error: " << e.what() << endl;
        }
        //
        // fast5 file is closed by its destructor at the end of this scope
        //
    }
    assert(hdf5_tools::File::get_object_count() == 0);
}
