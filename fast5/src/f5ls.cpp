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
            if (f.have_channel_id_params())
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
            if (f.have_tracking_id_params())
            {
                auto tracking_id_params = f.get_tracking_id_params();
                print_map(cout, tracking_id_params, "tracking_id/");
            }
            //
            // inspect sequences params
            //
            if (f.have_sequences_params())
            {
                auto sequences_params = f.get_sequences_params();
                print_map(cout, sequences_params, "sequences/");
            }
            //
            // inspect raw samples
            //
            if (f.have_raw_samples())
            {
                auto rs_params = f.get_raw_samples_params();
                auto rs = f.get_raw_samples();
                cout << "raw_samples/read_id=" << rs_params.read_id << endl
                     << "raw_samples/read_number=" << rs_params.read_number << endl
                     << "raw_samples/start_mux=" << rs_params.start_mux << endl
                     << "raw_samples/start_time=" << rs_params.start_time << endl
                     << "raw_samples/duration=" << rs_params.duration << endl
                     << "raw_samples/size=" << rs.size() << endl;
                const auto& e = rs.front();
                cout << "  (" << e << ")" << endl;
            }
            //
            // inspect eventdetection events
            //
            cout << "eventdetection_group_list=";
            print_vector(cout, f.get_eventdetection_group_list(), ",");
            cout << endl;
            if (f.have_eventdetection_events())
            {
                auto ed_params = f.get_eventdetection_params();
                print_map(cout, ed_params, "eventdetection/");
                auto ed_ev_params = f.get_eventdetection_events_params();
                auto ed_ev = f.get_eventdetection_events();
                cout << "eventdetection/events/abasic_found=" << ed_ev_params.abasic_found << endl
                     << "eventdetection/events/duration=" << ed_ev_params.duration << endl
                     << "eventdetection/events/median_before=" << ed_ev_params.median_before << endl
                     << "eventdetection/events/read_id=" << ed_ev_params.read_id << endl
                     << "eventdetection/events/read_number=" << ed_ev_params.read_number << endl
                     << "eventdetection/events/scaling_used=" << ed_ev_params.scaling_used << endl
                     << "eventdetection/events/start_mux=" << ed_ev_params.start_mux << endl
                     << "eventdetection/events/start_time=" << ed_ev_params.start_time << endl
                     << "eventdetection/events/size=" << ed_ev.size() << endl;
                const auto& e = ed_ev.front();
                cout << "  (mean=" << e.mean
                     << ", stdv=" << e.stdv
                     << ", start=" << e.start
                     << ", length=" << e.length << ")" << endl;
            } // if have_eventdetection_events
            //
            // inspect basecall groups
            //
            for (unsigned st = 0; st < 3; ++st)
            {
                cout << "basecall(" << st << ")/group_list=";
                print_vector(cout, f.get_basecall_strand_group_list(st), ",");
                cout << endl;
                // basecall sequence
                if (f.have_basecall_seq(st))
                {
                    cout << "basecall(" << st << ")/seq_size=" << f.get_basecall_seq(st).size() << endl;
                }
                // basecall model
                if (f.have_basecall_model(st))
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
                    const auto& e = m.front();
                    cout << "  (kmer=" << e.get_kmer()
                         << ", level_mean=" << e.level_mean
                         << ", level_stdv=" << e.level_stdv << ")" << endl;
                }
                // basecall events
                if (f.have_basecall_events(st))
                {
                    auto ev = f.get_basecall_events(st);
                    cout << "basecall(" << st << ")/events/size=" << ev.size() << endl;
                    const auto& e = ev.front();
                    cout << "  (mean=" << e.mean
                         << ", stdv=" << e.stdv
                         << ", start=" << e.start
                         << ", length=" << e.length
                         << ", model_state=" << e.get_model_state()
                         << ", p_model_state=" << e.p_model_state
                         << ", move=" << e.move << ")" << endl;
                }
                // basecall alignment
                if (st == 2 and f.have_basecall_alignment())
                {
                    auto al = f.get_basecall_alignment();
                    cout << "basecall(2)/alignment/size=" << al.size() << endl;
                    const auto& e = al.front();
                    cout << "  (template_index=" << e.template_index
                         << ", complement_index=" << e.complement_index
                         << ", kmer=" << e.get_kmer() << ")" << endl;
                }
            } // for st
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
