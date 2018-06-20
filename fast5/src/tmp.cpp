//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#include <cassert>
#include <exception>
#include <functional>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <deque>

#include <hdf5.h>

using namespace std;

template < typename T, typename U >
std::size_t offset_of(U T::* mem_ptr)
{
    return reinterpret_cast< std::size_t >(&(((T*)0)->*mem_ptr));
}

struct A
{
    int val_1;
    unsigned val_2;
    float val_3;
    int val_4;
    string val_5;
};

int main(int argc, char * argv[])
{
    if (argc != 2)
    {
        cerr << "use: " << argv[0] << " <file>" << endl;
        exit(EXIT_FAILURE);
    }
    //
    // create file, fail if existing
    //
    auto file_id = H5Fcreate(argv[1], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id > 0);
    auto scalar_space_id = H5Screate(H5S_SCALAR);
    assert(scalar_space_id > 0);
    auto lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl_id, 1);

    //
    // write numeric scalar attribute
    //
    // create group
    auto grp_id = H5Gcreate2(file_id, "/Group_1/Subgroup_1_1", lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
    assert(grp_id > 0);
    auto attr1_id = H5Acreate2(grp_id, "Attribute_1_1_1", H5T_NATIVE_INT, scalar_space_id,
                               H5P_DEFAULT, H5P_DEFAULT);
    assert(attr1_id > 0);
    int i = 42;
    auto status = H5Awrite(attr1_id, H5T_NATIVE_INT, &i);
    assert(status >= 0);
    H5Gclose(grp_id);
    H5Aclose(attr1_id);

    //
    // write numeric vector dataset
    //
    {
        vector< float > v = { 1.0, 2.0, 3.0 };
        hsize_t v_size = v.size();
        auto v_space_id = H5Screate_simple(1, &v_size, nullptr);
        auto ds1_id = H5Dcreate2(file_id, "/Group_2/Subgroup_2_1/Dataset_2_1_1", H5T_NATIVE_FLOAT, v_space_id,
                                 lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(ds1_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
        assert(status >= 0);
        H5Dclose(ds1_id);
        H5Sclose(v_space_id);
    }

    //
    // write compound scalar dataset
    //
    {
        //A a{ 1, 2, 3.2, 4, true };
        A a{ 1, 2, 3.2, 4, "xoxo" };
        auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(A));
        vector< hid_t > a_stype_id;
        status = H5Tinsert(a_type_id, "val_2", offset_of(&A::val_2), H5T_NATIVE_UINT);
        assert(status >= 0);
        status = H5Tinsert(a_type_id, "val_3", offset_of(&A::val_3), H5T_NATIVE_FLOAT);
        assert(status >= 0);
        status = H5Tinsert(a_type_id, "val_1", offset_of(&A::val_1), H5T_NATIVE_INT);
        assert(status >= 0);
        auto ds2_id = H5Dcreate2(file_id, "/Group_2/Subgroup_2_1/Dataset_2_1_2", a_type_id, scalar_space_id,
                                 lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(ds2_id, a_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a);
        assert(status >= 0);
        H5Dclose(ds2_id);
        H5Tclose(a_type_id);
    }

    //
    // write compound scalar dataset
    //
    {
        //vector< A > a_v{{ 1, 2, 3.1, 4, true }, { 11, 12, 13.1, 14, false }, { 21, 22, 23.1, 24, true }};
        vector< A > a_v{{ 1, 2, 3.1, 4, "xoxo" }, { 11, 12, 13.1, 14, "xexe" }, { 21, 22, 23.1, 24, "xixi" }};
        auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(A));
        status = H5Tinsert(a_type_id, "val_1", offset_of(&A::val_1), H5T_NATIVE_INT);
        assert(status >= 0);
        //status = H5Tinsert(a_type_id, "val_2", offset_of(&A::val_2), H5T_NATIVE_UINT);
        //assert(status >= 0);
        status = H5Tinsert(a_type_id, "val_3", offset_of(&A::val_3), H5T_NATIVE_FLOAT);
        assert(status >= 0);
        hsize_t a_v_size = a_v.size();
        auto a_v_space_id = H5Screate_simple(1, &a_v_size, nullptr);
        auto ds3_id = H5Dcreate2(file_id, "/Group_2/Subgroup_2_1/Dataset_2_1_3", a_type_id, a_v_space_id,
                                 lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(ds3_id, a_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_v.data());
        assert(status >= 0);
        H5Dclose(ds3_id);
        H5Sclose(a_v_space_id);
        H5Tclose(a_type_id);
    }

    //
    // write compound scalar dataset in 2 steps
    //
    {
        //vector< A > a_v{{ 1, 2, 3.1, 4, true }, { 11, 12, 13.1, 14, false }, { 21, 22, 23.1, 24, true }};
        vector< A > a_v{{ 100, 2, 3.1, 4, "xoxo" }, { 111, 12, 13.1, 14, "xexe" }, { 121, 22, 23.1, 24, "xixi" }};
        hid_t ds4_id;
        // create dataset
        {
            auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(A));
            status = H5Tinsert(a_type_id, "val_1", offset_of(&A::val_1), H5T_NATIVE_INT);
            assert(status >= 0);
            status = H5Tinsert(a_type_id, "val_2", offset_of(&A::val_2), H5T_NATIVE_UINT);
            assert(status >= 0);
            status = H5Tinsert(a_type_id, "val_3", offset_of(&A::val_3), H5T_NATIVE_FLOAT);
            assert(status >= 0);
            hid_t val_5_type_id = H5Tcopy(H5T_C_S1);
            status = H5Tset_size(val_5_type_id, H5T_VARIABLE);
            assert(status >= 0);
            status = H5Tinsert(a_type_id, "val_5", offset_of(&A::val_5), val_5_type_id);
            assert(status >= 0);
            hsize_t a_v_size = a_v.size();
            auto a_v_space_id = H5Screate_simple(1, &a_v_size, nullptr);
            ds4_id = H5Dcreate2(file_id, "/Group_2/Subgroup_2_1/Dataset_2_1_4", a_type_id, a_v_space_id,
                                lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(a_v_space_id);
            H5Tclose(val_5_type_id);
            H5Tclose(a_type_id);
        }
        // write val_1
        {
            auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(A));
            status = H5Tinsert(a_type_id, "val_1", offset_of(&A::val_1), H5T_NATIVE_INT);
            assert(status >= 0);
            status = H5Dwrite(ds4_id, a_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_v.data());
            assert(status >= 0);
            H5Tclose(a_type_id);
        }
        // write val_2
        {
            auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(A));
            status = H5Tinsert(a_type_id, "val_2", offset_of(&A::val_2), H5T_NATIVE_UINT);
            assert(status >= 0);
            status = H5Dwrite(ds4_id, a_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_v.data());
            assert(status >= 0);
            H5Tclose(a_type_id);
        }
        // write val_3
        {
            auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(A));
            status = H5Tinsert(a_type_id, "val_3", offset_of(&A::val_3), H5T_NATIVE_FLOAT);
            assert(status >= 0);
            //status = H5Dwrite(ds4_id, a_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_v.data());
            assert(status >= 0);
            H5Tclose(a_type_id);
        }
        // write val_5
        {
            auto a_type_id = H5Tcreate(H5T_COMPOUND, sizeof(const char *));
            hid_t val_5_type_id = H5Tcopy(H5T_C_S1);
            status = H5Tset_size(val_5_type_id, H5T_VARIABLE);
            assert(status >= 0);
            status = H5Tinsert(a_type_id, "val_5", 0, val_5_type_id);
            assert(status >= 0);
            H5Tclose(val_5_type_id);
            vector< const char * > charptr_buff(a_v.size());
            for (size_t i = 0; i < a_v.size(); ++i)
            {
                charptr_buff[i] = a_v[i].val_5.data();
            }
            status = H5Dwrite(ds4_id, a_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, charptr_buff.data());
            assert(status >= 0);
            H5Tclose(a_type_id);
        }
        H5Dclose(ds4_id);
    }

    //
    // clean up
    //
    H5Sclose(scalar_space_id);
    H5Pclose(lcpl_id);
    H5Fclose(file_id);
}
