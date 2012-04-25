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

#define BOOST_TEST_MODULE accumulating_value
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>

#include <halmd/utility/program_options/program_options.hpp>
#include <test/tools/ctest.hpp>
#include <test/unit/utility/program_options/predicates.hpp>

using namespace boost;
using namespace std;

namespace po = halmd::po;

/**
 * test empty accumulating value
 */
BOOST_AUTO_TEST_CASE( accumulating_value_empty )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK( vm["verbose"].empty() );
}

/**
 * test accumulating value with 1 occurrence
 */
BOOST_AUTO_TEST_CASE( accumulating_value_1 )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "positional", "--verbose"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 1 );
}

/**
 * test accumulating value with 3 occurrences
 */
BOOST_AUTO_TEST_CASE( accumulating_value_3 )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "-v", "positional", "--verbose", "-v"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 3 );
}

/**
 * test default accumulating value
 */
BOOST_AUTO_TEST_CASE( accumulating_value_defaulted )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>()->default_value(2), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), 2 );
    BOOST_CHECK( vm["verbose"].defaulted() );
}

/**
 * test accumulating value with defaulted value and many occurrences
 */
BOOST_AUTO_TEST_CASE( accumulating_value_defaulted_many )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>()->default_value(3), "")
        ;
    vector<char const*> args;
    args.push_back(""); //< argv[0]
    for (size_t i = 1; i <= 42; ++i) {
        args.push_back("-v");
        args.push_back("positional");
        args.push_back("--verbose");
        args.push_back("-v");
    }
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), (42 * 3) + 3 );
    BOOST_CHECK( !vm["verbose"].defaulted() );
}

/**
 * test conflicting accumulating value with defaulted value and many occurrences
 */
BOOST_AUTO_TEST_CASE( conflicting_accumulating_value_defaulted_many )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>()->default_value(3), "")
        ("verbosity", po::accum_value<int>()->conflicts("verbose")->default_value(4), "")
        ;
    vector<char const*> args;
    args.push_back(""); //< argv[0]
    for (size_t i = 1; i <= 42; ++i) {
        args.push_back("-v");
        args.push_back("positional");
        args.push_back("--verbose");
        args.push_back("-v");
    }
    args.push_back("--verbosity");
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::conflicting_option>(parsed, vm, "verbosity") );
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), (42 * 3) + 3 );
    BOOST_CHECK( !vm["verbose"].defaulted() );
    BOOST_CHECK_EQUAL( vm["verbosity"].as<int>(), 4 + 1 );
    BOOST_CHECK( !vm["verbosity"].defaulted() );
}

/**
 * test dependent accumulating value with defaulted value and many occurrences
 */
BOOST_AUTO_TEST_CASE( dependent_accumulating_value_defaulted_many )
{
    po::options_description desc;
    desc.add_options()
        ("verbose,v", po::accum_value<int>()->depends("verbosity")->default_value(3), "")
        ("verbosity", po::accum_value<int>()->default_value(4)->conflicts("verbose"), "")
        ;
    vector<char const*> args;
    args.push_back(""); //< argv[0]
    for (size_t i = 1; i <= 42; ++i) {
        args.push_back("-v");
        args.push_back("positional");
        args.push_back("--verbose");
        args.push_back("-v");
    }
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::dependent_option>(parsed, vm, "verbose") );
    BOOST_CHECK_NO_THROW( po::notify(vm) );
    BOOST_CHECK_EQUAL( vm["verbose"].as<int>(), (42 * 3) + 3 );
    BOOST_CHECK( !vm["verbose"].defaulted() );
    BOOST_CHECK_EQUAL( vm["verbosity"].as<int>(), 4 );
    BOOST_CHECK( vm["verbosity"].defaulted() );
}
