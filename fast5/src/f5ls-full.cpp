//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#include <cassert>
#include <iostream>
#include <string>

#include "fast5.hpp"

using namespace std;

template < typename T >
void print_vector(ostream& os, const vector< T >& v, const string& delim)
{
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        if (it != v.begin()) os << delim;
        os << *it;
    }
}
template < typename U, typename V >
void print_map(ostream& os, const map< U, V >& m, const string& prefix)
{
    for (const auto& p : m)
    {
        os << prefix << p.first << "=" << p.second << endl;
    }
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "use: " << argv[0] << " <fast5_file>" << endl;
        return EXIT_FAILURE;
    }
    string file_name(argv[1]);
    //
    // open the FAST5 file for reading
    //
    if (not fast5::File::is_valid_file(file_name))
    {
        cout << "not a fast5 file [" << file_name << "]" << endl;
        return EXIT_SUCCESS;
    }
    {
        fast5::File f;
        //
        // All fast5 operations are performed inside a try-catch block. This should
        // resist various hdf5 errors without leaking memory.
        //
        try
        {
            //
            // open file
            //
            f.open(file_name);
            assert(f.is_open());
            //
            // extract version information for the ONT software used to generate this dataset
            //
            cout << "file_version=" << f.file_version() << endl;
            //
            // inspect channel_id params
            //
            bool have_channel_id_params = f.have_channel_id_params();
            cout << "have_channel_id_params=" << have_channel_id_params << endl;
            if (have_channel_id_params)
            {
                auto channel_id_params = f.get_channel_id_params();
                cout << "channel_id/channel_number=" << channel_id_params.channel_number << endl
                     << "channel_id/digitisation=" << channel_id_params.digitisation << endl
                     << "channel_id/offset=" << channel_id_params.offset << endl
                     << "channel_id/range=" << channel_id_params.range << endl
                     << "channel_id/sampling_rate=" << channel_id_params.sampling_rate << endl;
            }
            //
            // inspect tracking_id params
            //
            bool have_tracking_id_params = f.have_tracking_id_params();
            cout << "have_tracking_id_params=" << have_tracking_id_params << endl;
            if (have_tracking_id_params)
            {
                auto tracking_id_params = f.get_tracking_id_params();
                print_map(cout, tracking_id_params, "tracking_id/");
            }
            //
            // inspect sequences params
            //
            bool have_sequences_params = f.have_sequences_params();
            cout << "have_sequences_params=" << have_sequences_params << endl;
            if (have_sequences_params)
            {
                auto sequences_params = f.get_sequences_params();
                print_map(cout, sequences_params, "sequences/");
            }
            //
            // inspect raw samples
            //
            bool have_raw_samples = f.have_raw_samples();
            cout << "have_raw_samples=" << have_raw_samples << endl;
            if (have_raw_samples)
            {
                auto rs_rn_list = f.get_raw_samples_read_name_list();
                cout << "raw_samples_read_name_list=";
                print_vector(cout, rs_rn_list, ",");
                cout << endl;
                for (const auto& rn : rs_rn_list)
                {
                    auto rs_params = f.get_raw_samples_params();
                    auto rs = f.get_raw_samples();
                    cout << "raw_samples/" << rn << "/read_id=" << rs_params.read_id << endl
                         << "raw_samples/" << rn << "/read_number=" << rs_params.read_number << endl
                         << "raw_samples/" << rn << "/start_mux=" << rs_params.start_mux << endl
                         << "raw_samples/" << rn << "/start_time=" << rs_params.start_time << endl
                         << "raw_samples/" << rn << "/duration=" << rs_params.duration << endl
                         << "raw_samples/" << rn << "/size=" << rs.size() << endl;
                    const auto& e = rs.front();
                    cout << "  (" << e << ")" << endl;
                }
            }
            //
            // inspect eventdetection groups
            //
            bool have_eventdetection_events = f.have_eventdetection_events();
            cout << "have_eventdetection_events=" << have_eventdetection_events << endl;
            bool have_eventdetection_group = f.have_eventdetection_group();
            cout << "have_eventdetection_group=" << have_eventdetection_group << endl;
            if (have_eventdetection_group)
            {
                auto ed_gr_list = f.get_eventdetection_group_list();
                cout << "eventdetection_group_list=";
                print_vector(cout, ed_gr_list, ",");
                cout << endl;
                for (const auto& ed_gr : ed_gr_list)
                {
                    auto ed_params = f.get_eventdetection_params(ed_gr);
                    print_map(cout, ed_params, "eventdetection/");
                    auto rn_list = f.get_eventdetection_read_name_list(ed_gr);
                    cout << "eventdetection/" << ed_gr << "/read_name_list=";
                    print_vector(cout, rn_list, ",");
                    cout << endl;
                    have_eventdetection_events = f.have_eventdetection_events(ed_gr);
                    cout << "eventdetection/" << ed_gr << "/have_eventdetection_events=" << have_eventdetection_events << endl;
                    for (const auto& rn : rn_list)
                    {
                        std::ostringstream tmp;
                        tmp << "eventdetection/" << ed_gr << "/" << rn;
                        auto ed_ev_params = f.get_eventdetection_events_params(ed_gr, rn);
                        auto ed_ev = f.get_eventdetection_events(ed_gr, rn);
                        cout << tmp.str() << "/abasic_found=" << ed_ev_params.abasic_found << endl
                             << tmp.str() << "/duration=" << ed_ev_params.duration << endl
                             << tmp.str() << "/median_before=" << ed_ev_params.median_before << endl
                             << tmp.str() << "/read_id=" << ed_ev_params.read_id << endl
                             << tmp.str() << "/read_number=" << ed_ev_params.read_number << endl
                             << tmp.str() << "/scaling_used=" << ed_ev_params.scaling_used << endl
                             << tmp.str() << "/start_mux=" << ed_ev_params.start_mux << endl
                             << tmp.str() << "/start_time=" << ed_ev_params.start_time << endl
                             << tmp.str() << "/size=" << ed_ev.size() << endl;
                        for (const auto& e : ed_ev)
                        {
                            cout << "  (mean=" << e.mean
                                 << ", stdv=" << e.stdv
                                 << ", start=" << e.start
                                 << ", length=" << e.length << ")" << endl;
                            break;
                        }
                    } // for rn : rn_list
                } // for ed_gr : ed_gr_list
            } // if have_eventdetection_groups
            //
            // inspect basecall groups
            //
            bool have_basecall_group = f.have_basecall_group();
            cout << "have_basecall_group=" << have_basecall_group << endl;
            if (have_basecall_group)
            {
                auto bc_gr_list = f.get_basecall_group_list();
                cout << "basecall_group_list=";
                print_vector(cout, bc_gr_list, ",");
                cout << endl;
                for (unsigned st = 0; st < 3; ++st)
                {
                    auto bc_st_gr_list = f.get_basecall_strand_group_list(st);
                    cout << "basecall_strand_group_list(" << st << ")=";
                    print_vector(cout, bc_st_gr_list, ",");
                    cout << endl;
                }
                for (const auto& bc_gr : bc_gr_list)
                {
                    // dump basecall params
                    auto bc_params = f.get_basecall_params(bc_gr);
                    std::ostringstream tmp;
                    tmp << "basecall/" << bc_gr << "/";
                    print_map(cout, bc_params, tmp.str());
                    // check if basecall log exists
                    cout << "basecall/" << bc_gr << "/have_log=" << f.have_basecall_log(bc_gr) << endl;
                    // check if eventdetection link exists
                    auto bc_ed_gr = f.get_basecall_eventdetection_group(bc_gr);
                    cout << "basecall/" << bc_gr << "/eventdetection_group=" << bc_ed_gr << endl;
                }
                for (unsigned st = 0; st < 3; ++st)
                {
                    auto gr_l = f.get_basecall_strand_group_list(st);
                    bool have_seq = f.have_basecall_seq(st);
                    cout << "basecall(" << st << ")/have_seq=" << have_seq << endl;
                    if (have_seq)
                    {
                        cout << "basecall(" << st << ")/seq=" << f.get_basecall_seq(st).substr(0, 10) << "..." << endl;
                    }
                    bool have_model = f.have_basecall_model(st);
                    cout << "basecall(" << st << ")/have_model=" << have_model << endl;
                    if (have_model)
                    {
                        cout << "basecall(" << st << ")/model_file=" << f.get_basecall_model_file(st) << endl;
                        auto m_params = f.get_basecall_model_params(st);
                        auto m = f.get_basecall_model(st);
                        cout << "basecall(" << st << ")/model/scale=" << m_params.scale << endl
                             << "basecall(" << st << ")/model/shift=" << m_params.shift << endl
                             << "basecall(" << st << ")/model/drift=" << m_params.drift << endl
                             << "basecall(" << st << ")/model/var=" << m_params.var << endl
                             << "basecall(" << st << ")/model/scale_sd=" << m_params.scale_sd << endl
                             << "basecall(" << st << ")/model/var_sd=" << m_params.var_sd << endl
                             << "basecall(" << st << ")/model/size=" << m.size() << endl;
                        for (const auto& e : m)
                        {
                            cout << "  (kmer=" << e.get_kmer()
                                 << ", level_mean=" << e.level_mean
                                 << ", level_stdv=" << e.level_stdv << ")" << endl;
                            break;
                        }
                    }
                    bool have_events = f.have_basecall_events(st);
                    cout << "basecall(" << st << ")/have_events=" << have_events << endl;
                    if (have_events)
                    {
                        auto ev = f.get_basecall_events(st);
                        cout << "basecall(" << st << ")/events/size=" << ev.size() << endl;
                        for (const auto& e : ev)
                        {
                            cout << "  (mean=" << e.mean
                                 << ", stdv=" << e.stdv
                                 << ", start=" << e.start
                                 << ", length=" << e.length
                                 << ", model_state=" << e.get_model_state()
                                 << ", p_model_state=" << e.p_model_state
                                 << ", move=" << e.move << ")" << endl;
                            break;
                        }
                    }
                    if (st == 2)
                    {
                        bool have_alignment = f.have_basecall_alignment();
                        cout << "basecall(2)/have_alignment=" << have_alignment << endl;
                        if (have_alignment)
                        {
                            auto al = f.get_basecall_alignment();
                            cout << "basecall(2)/alignment/size=" << al.size() << endl;
                            for (const auto& e : al)
                            {
                                cout << "  (template_index=" << e.template_index
                                     << ", complement_index=" << e.complement_index
                                     << ", kmer=" << e.get_kmer() << ")" << endl;
                                break;
                            }
                        }
                    }
                }
            } // have_basecall_groups
        }
        catch (hdf5_tools::Exception& e)
        {
            cout << "hdf5 error: " << e.what() << endl;
        }
        //
        // fast5 file is closed by its destructor at the end of this scope
        //
    }
    assert(fast5::File::get_object_count() == 0);
}
