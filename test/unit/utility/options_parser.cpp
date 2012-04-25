/*
 * Copyright © 2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE options_parser
#include <boost/test/unit_test.hpp>

#include <boost/version.hpp>
#include <sstream>

#include <halmd/utility/options_parser.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/lua.hpp>

using namespace boost;
using namespace std;

namespace po = halmd::po;

// macros do not accept comma-separate template parameters
typedef multi_array<unsigned int, 1> uint_array;
typedef multi_array<double, 1> double_array;

BOOST_AUTO_TEST_SUITE( cplusplus )

BOOST_AUTO_TEST_SUITE( command_line_parser )

/**
 * test empty arguments
 */
BOOST_AUTO_TEST_CASE( empty_args )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("particles", po::value<uint_array>(), "")
        ;
    array<char const*, 1> args = {{ "halmd" //< argv[0]
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 0LU );
}

/**
 * test empty arguments with defaulted value
 */
BOOST_AUTO_TEST_CASE( empty_args_defaulted )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("density", po::value<double>()->default_value(0.75), "")
        ;
    array<char const*, 1> args = {{ "halmd" //< argv[0]
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.75 );
}

/**
 * test empty options description
 */
BOOST_AUTO_TEST_CASE( empty_options_description )
{
    halmd::options_parser parser;
    po::options_description desc;
    parser.add(desc, "core");
    array<char const*, 2> args = {{ "halmd" //< argv[0]
      , "core"
    }};
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm), po::too_many_positional_options_error );
}

/**
 * test single occurrence of module namespace
 */
BOOST_AUTO_TEST_CASE( single_occurrence )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("particles", po::value<uint_array>(), "")
        ("density", po::value<double>(), "")
        ;
    array<char const*, 6> args = {{ "" //< argv[0]
      , "box", "--particles", "10000" , "--density", "0.8"
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
}

/**
 * test multiple occurrences of module namespace
 */
BOOST_AUTO_TEST_CASE( multiple_occurrences )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("particles", po::value<uint_array>(), "")
        ("density", po::value<double>(), "")
        ;
    array<char const*, 7> args = {{ "" //< argv[0]
      , "box", "--particles", "10000"
      , "box", "--density", "0.8"
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
}

/**
 * test module-independent options
 */
BOOST_AUTO_TEST_CASE( module_independent_options )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>()->default_value(42), "")
        ("output", po::value<string>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "-vvv", "-v", "--output", "/dev/null"
    }};
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 2LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), (42 + 4) );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "/dev/null" );
}

/**
 * test module-independent and module options
 */
BOOST_AUTO_TEST_CASE( module_independent_and_module_options )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    po::options_description desc;
    desc.add_options()
        ("particles", po::value<uint_array>(), "")
        ("density", po::value<double>(), "")
        ;
    array<char const*, 9> args = {{ "" //< argv[0]
      , "-vvv", "-v"
      , "box", "--particles", "10000"
      , "box", "--density", "0.8"
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 3LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 4 );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "halmd" );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
}

/**
 * test multiple occurrences of option
 */
BOOST_AUTO_TEST_CASE( multiple_option_occurrences )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("dimension", po::value<unsigned int>(), "")
        ;
    array<char const*, 4> args = {{ "" //< argv[0]
      , "box", "--dimension=2", "--dimension=3"
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm), po::multiple_occurrences);
}

/**
 * test multiple occurrences of namespace and option
 */
BOOST_AUTO_TEST_CASE( multiple_namespace_and_option_occurrences )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("integrator", po::value<string>(), "")
            ;
        parser.add(desc, "core");
    }
    array<char const*, 7> args = {{ "" //< argv[0]
      , "box", "--dimension=2"
      , "core", "--integrator=verlet_nvt_andersen"
      , "box", "--dimension=3"
    }};
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm), po::multiple_occurrences);
}

/**
 * test invalid module namespace
 */
BOOST_AUTO_TEST_CASE( invalid_namespace )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("dimension", po::value<unsigned int>(), "")
        ;
    array<char const*, 4> args = {{ "" //< argv[0]
      , "box", "--dimension=2"
      , "core"
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm), po::too_many_positional_options_error );
}

/**
 * test invalid module namespace with option
 */
BOOST_AUTO_TEST_CASE( invalid_namespace_with_option )
{
    halmd::options_parser parser;
    po::options_description desc;
    desc.add_options()
        ("dimension", po::value<unsigned int>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "box", "--dimension=2"
      , "core", "--integrator=verlet_nvt_andersen"
    }};
    parser.add(desc, "box");
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm), po::unknown_option );
}

/**
 * test multiple modules
 */
BOOST_AUTO_TEST_CASE( multiple_modules )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("integrator", po::value<string>(), "")
            ("force", po::value<string>(), "")
            ;
        parser.add(desc, "core");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("smooth", po::value<double>()->implicit_value(0.001), "")
            ;
        parser.add(desc, "morse");
    }
    array<char const*, 15> args = {{ "" //< argv[0]
      , "-vvv"
      , "box", "--dimension=3"
      , "core", "--integrator=verlet_nvt_andersen"
      , "box", "--particles", "10000" , "--density", "0.8"
      , "morse", "--smooth"
      , "core", "--force=morse"
    }};
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    BOOST_CHECK_EQUAL( vm.size(), 5LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 3 );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "halmd" );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["dimension"].as<unsigned int>(), 3U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["smooth"].as<double>(), 0.001 );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["integrator"].as<string>(), "verlet_nvt_andersen" );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["force"].as<string>(), "morse" );
}

/**
 * test multiple modules using std::vector for arguments
 */
BOOST_AUTO_TEST_CASE( multiple_modules_with_vector_args )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("integrator", po::value<string>(), "")
            ("force", po::value<string>(), "")
            ;
        parser.add(desc, "core");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("smooth", po::value<double>()->implicit_value(0.001), "")
            ;
        parser.add(desc, "morse");
    }
    vector<string> args;
    args.push_back("-vvv");
    args.push_back("box");
    args.push_back("--dimension=3");
    args.push_back("core");
    args.push_back("--integrator=verlet_nvt_andersen");
    args.push_back("box");
    args.push_back("--particles");
    args.push_back("10000" );
    args.push_back("--density");
    args.push_back("0.8");
    args.push_back("morse");
    args.push_back("--smooth");
    args.push_back("core");
    args.push_back("--force=morse");
    po::variables_map vm;
    parser.parse_command_line(args, vm);
    BOOST_CHECK_EQUAL( vm.size(), 5LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 3 );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "halmd" );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["dimension"].as<unsigned int>(), 3U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["smooth"].as<double>(), 0.001 );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["integrator"].as<string>(), "verlet_nvt_andersen" );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["force"].as<string>(), "morse" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( config_file_parser )

/**
 * test empty arguments
 */
BOOST_AUTO_TEST_CASE( empty_args )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 0LU );
}

/**
 * test empty arguments with defaulted value
 */
BOOST_AUTO_TEST_CASE( empty_args_defaulted )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("density", po::value<double>()->default_value(0.75), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.75 );
}

/**
 * test unknown option
 */
BOOST_AUTO_TEST_CASE( unknown_option )
{
    halmd::options_parser parser;
    parser.add_options()
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    stringstream is;
    is << "file  = /dev/null" << endl
       ;
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_config_file(is, vm), po::unknown_option);
}

/**
 * test unknown module option
 */
BOOST_AUTO_TEST_CASE( unknown_module_option )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "[box]"               << endl
       << "Particles = 10000"   << endl
       ;
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_config_file(is, vm), po::unknown_option);
}

boost::test_tools::predicate_result
parse_config_file_throws_error(
    halmd::options_parser& parser
  , std::istream& is, boost::program_options::variables_map& vm
  , std::string const& error)
{
    try {
        parser.parse_config_file(is, vm);
    }
    catch (boost::program_options::error const& e) {
        if (e.what() != error) {
            boost::test_tools::predicate_result res(false);
            res.message()
                << "Exception has wrong error [" << e.what()
                << " != " << error << "]"
                ;
            return res;
        }
        return true;
    }
    boost::test_tools::predicate_result res(false);
    res.message() << "No exception thrown";
    return res;
}

/**
 * test empty options description
 */
BOOST_AUTO_TEST_CASE( empty_options_description )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        parser.add(desc, "core");
    }
    stringstream is;
    is << "[core]"              << endl
       << "force = morse"       << endl;
    po::variables_map vm;
    BOOST_CHECK( parse_config_file_throws_error(parser, is, vm, "unknown module " "core") );
}

/**
 * test single occurrence of module namespace
 */
BOOST_AUTO_TEST_CASE( single_occurrence )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "[box]"               << endl
       << "particles = 10000"   << endl
       << "density = 0.8"       << endl
       ;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
}

/**
 * test multiple occurrences of module namespace
 */
BOOST_AUTO_TEST_CASE( multiple_occurrences )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "[box]"               << endl
       << "particles = 10000"   << endl
       << "[box]"               << endl
       << "density = 0.8"       << endl
       ;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
}

/**
 * test module-independent options
 */
BOOST_AUTO_TEST_CASE( module_independent_options )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>()->default_value(42), "")
        ("output", po::value<string>(), "")
        ;
    stringstream is;
    is << "verbose = true"      << endl
       << "verbose = true"      << endl
       << "verbose = true"      << endl
       << "output  = /dev/null" << endl
       << "verbose = true"      << endl
       ;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 2LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), (42 + 4) );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "/dev/null" );
}

/**
 * test module-independent and module options
 */
BOOST_AUTO_TEST_CASE( module_independent_and_module_options )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "verbose = true"      << endl
       << "verbose = true"      << endl
       << "verbose = true"      << endl
       << "verbose = true"      << endl
       << "[box]"               << endl
       << "particles = 10000"   << endl
       << "[box]"               << endl
       << "density = 0.8"       << endl
       ;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 3LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 4 );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "halmd" );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
}

/**
 * test multiple occurrences of option
 */
BOOST_AUTO_TEST_CASE( multiple_option_occurrences )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "[box]"               << endl
       << "dimension=2"         << endl
       << "dimension=3"         << endl
       ;
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_config_file(is, vm), po::multiple_occurrences);
}

/**
 * test multiple occurrences of namespace and option
 */
BOOST_AUTO_TEST_CASE( multiple_namespace_and_option_occurrences )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("integrator", po::value<string>(), "")
            ;
        parser.add(desc, "core");
    }
    stringstream is;
    is << "[box]"                               << endl
       << "dimension=2"                         << endl
       << "[core]"                              << endl
       << "integrator=verlet_nvt_andersen"      << endl
       << "[box]"                               << endl
       << "dimension=3"                         << endl
       ;
    po::variables_map vm;
    BOOST_CHECK_THROW( parser.parse_config_file(is, vm), po::multiple_occurrences);
}

/**
 * test invalid module namespace
 */
BOOST_AUTO_TEST_CASE( invalid_namespace )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "[box]"               << endl
       << "dimension=2"         << endl
       << "[core]"              << endl
       ;
    po::variables_map vm;
    BOOST_CHECK_NO_THROW( parser.parse_config_file(is, vm) );
}

/**
 * test invalid module namespace with option
 */
BOOST_AUTO_TEST_CASE( invalid_namespace_with_option )
{
    halmd::options_parser parser;
    {
        po::options_description desc;
        desc.add_options()
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    stringstream is;
    is << "[box]"                               << endl
       << "dimension=2"                         << endl
       << "[core]"                              << endl
       << "integrator=verlet_nvt_andersen"      << endl
       ;
    po::variables_map vm;
    BOOST_CHECK( parse_config_file_throws_error(parser, is, vm, "unknown module " "core") );
}

/**
 * test multiple modules
 */
BOOST_AUTO_TEST_CASE( multiple_modules )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("integrator", po::value<string>(), "")
            ("force", po::value<string>(), "")
            ;
        parser.add(desc, "core");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("epsilon", po::value<double_array>(), "")
            ;
        parser.add(desc, "morse");
    }
    stringstream is;
#if BOOST_VERSION >= 104200
    is << "verbose ="                           << endl
       << "verbose ="                           << endl
       << "verbose ="                           << endl
#else // https://svn.boost.org/trac/boost/ticket/1537
    is << "verbose = true"                      << endl
       << "verbose = true"                      << endl
       << "verbose = true"                      << endl
#endif
                                                << endl
       << "[box]"                               << endl
       << "dimension=3"                         << endl
       << "[core]"                              << endl
       << "integrator=verlet_nvt_andersen"      << endl
                                                << endl
       << "[box]"                               << endl
       << "particles = 10000"                   << endl
       << "density = 0.8"                       << endl
       << "[morse]"                             << endl
       << "epsilon = 1.0,0.88,0.8"              << endl
                                                << endl
       << "[core]"                              << endl
       << "force=morse"                         << endl
       ;
    po::variables_map vm;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 5LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 3 );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "halmd" );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["dimension"].as<unsigned int>(), 3U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 10000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>()[0], 1.0 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>()[1], 0.88 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>()[2], 0.8 );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["integrator"].as<string>(), "verlet_nvt_andersen" );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["force"].as<string>(), "morse" );
}

/**
 * test command line parser together with config file parser
 */
BOOST_AUTO_TEST_CASE( command_line_and_config_file )
{
    halmd::options_parser parser;
    parser.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ("output", po::value<string>()->default_value("halmd"), "")
        ;
    {
        po::options_description desc;
        desc.add_options()
            ("particles", po::value<uint_array>(), "")
            ("density", po::value<double>(), "")
            ("dimension", po::value<unsigned int>(), "")
            ;
        parser.add(desc, "box");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("integrator", po::value<string>(), "")
            ("force", po::value<string>(), "")
            ;
        parser.add(desc, "core");
    }
    {
        po::options_description desc;
        desc.add_options()
            ("epsilon", po::value<double_array>(), "")
            ("smooth", po::value<double>()->implicit_value(0.001), "")
            ;
        parser.add(desc, "morse");
    }
    array<char const*, 13> args = {{ "" //< argv[0]
      , "-v"
      , "box", "--dimension=2"
      , "core", "--integrator=verlet_nvt_andersen"
      , "box", "--particles", "20000"
      , "morse", "--smooth"
      , "core", "--force=morse"
    }};
    po::variables_map vm;
    parser.parse_command_line(args.size(), const_cast<char**>(&args.front()), vm);
    stringstream is;
#if BOOST_VERSION >= 104200
    is << "verbose ="                           << endl
       << "verbose ="                           << endl
       << "verbose ="                           << endl
#else // https://svn.boost.org/trac/boost/ticket/1537
    is << "verbose = true"                      << endl
       << "verbose = true"                      << endl
       << "verbose = true"                      << endl
#endif
                                                << endl
       << "[box]"                               << endl
       << "dimension=3"                         << endl
       << "[core]"                              << endl
       << "integrator=verlet_nvt_andersen"      << endl
                                                << endl
       << "[box]"                               << endl
       << "particles = 10000"                   << endl
       << "density = 0.8"                       << endl
       << "[morse]"                             << endl
       << "epsilon = 1.0,0.88,0.8"              << endl
                                                << endl
       << "[core]"                              << endl
       << "force=lennard_jones"                 << endl
       ;
    parser.parse_config_file(is, vm);
    BOOST_CHECK_EQUAL( vm.size(), 5LU );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 1 );
    BOOST_CHECK_EQUAL( vm["output"].as<string>(), "halmd" );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["dimension"].as<unsigned int>(), 2U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>().size(), 1LU );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["particles"].as<uint_array>()[0], 20000U );
    BOOST_CHECK_EQUAL( vm["box"].as<po::variables_map>()["density"].as<double>(), 0.8 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>()[0], 1.0 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>()[1], 0.88 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["epsilon"].as<double_array>()[2], 0.8 );
    BOOST_CHECK_EQUAL( vm["morse"].as<po::variables_map>()["smooth"].as<double>(), 0.001 );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["integrator"].as<string>(), "verlet_nvt_andersen" );
    BOOST_CHECK_EQUAL( vm["core"].as<po::variables_map>()["force"].as<string>(), "morse" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( lua )

// FIXME test halmd.modules.options using dummy Lua modules

BOOST_AUTO_TEST_SUITE_END()
