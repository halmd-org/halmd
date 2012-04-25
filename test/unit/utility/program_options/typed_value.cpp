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

#define BOOST_TEST_MODULE program_options
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>
#include <boost/version.hpp>
#include <stdint.h> // <cstdint> is C++0x

#include <halmd/utility/program_options/program_options.hpp>
#include <test/tools/ctest.hpp>
#include <test/unit/utility/program_options/predicates.hpp>

using namespace boost;
using namespace std;

namespace po = halmd::po;

/**
 * test conflicting option
 */
BOOST_AUTO_TEST_CASE( conflicting_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->conflicts("steps"), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296" //< MAX_UINT + 1LLU
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::conflicting_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test reversed conflicting option
 */
BOOST_AUTO_TEST_CASE( reversed_conflicting_option )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->conflicts("steps"), "")
        ("steps", po::value<uint64_t>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296"
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::conflicting_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test conflicting option with reversed arguments
 */
BOOST_AUTO_TEST_CASE( conflicting_option_reversed_args )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->conflicts("steps"), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--time", "123.4567"
      , "--steps", "4294967296"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::conflicting_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test reversed conflicting option with reversed arguments
 */
BOOST_AUTO_TEST_CASE( reversed_conflicting_option_reversed_args )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->conflicts("steps"), "")
        ("steps", po::value<uint64_t>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--time", "123.4567"
      , "--steps", "4294967296"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::conflicting_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test many conflicting option
 */
BOOST_AUTO_TEST_CASE( many_conflicting_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("seconds", po::value<double>(), "")
        ("time", po::value<double>()->conflicts("steps")->conflicts("seconds"), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--seconds", "123.4567e-15"
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::conflicting_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK( vm["steps"].empty() );
    BOOST_CHECK_EQUAL( vm["seconds"].as<double>(), 123.4567e-15 );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test non-conflicting option
 */
BOOST_AUTO_TEST_CASE( non_conflicting_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("seconds", po::value<double>(), "")
        ("time", po::value<double>()->conflicts("steps")->conflicts("seconds"), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK( vm["steps"].empty() );
    BOOST_CHECK( vm["seconds"].empty() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test conflicting option with defaulted values
 */
BOOST_AUTO_TEST_CASE( conflicting_option_defaulted )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>()->default_value(4294967296LLU), "")
        ("time", po::value<double>()->conflicts("steps")->default_value(123.4567), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["steps"].defaulted() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( vm["time"].defaulted() );
}

/**
 * test conflicting option with defaulted first value
 */
BOOST_AUTO_TEST_CASE( conflicting_option_defaulted_first )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>()->default_value(4294967296LLU), "")
        ("time", po::value<double>()->conflicts("steps"), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["steps"].defaulted() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( !vm["time"].defaulted() );
}

/**
 * test conflicting option with defaulted second value
 */
BOOST_AUTO_TEST_CASE( conflicting_option_defaulted_second )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->conflicts("steps")->default_value(123.4567), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--steps", "4294967296"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( !vm["steps"].defaulted() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( vm["time"].defaulted() );
}

/**
 * test non-existent conflicting option
 */
BOOST_AUTO_TEST_CASE( non_existent_conflicting_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->conflicts("step"), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296" //< MAX_UINT + 1LLU
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::unknown_option>(parsed, vm, "step") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test non-existent unused conflicting option
 */
BOOST_AUTO_TEST_CASE( non_existent_unused_conflicting_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("seconds", po::value<double>()->conflicts("step"), "")
        ("time", po::value<double>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296" //< MAX_UINT + 1LLU
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::unknown_option>(parsed, vm, "step") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["seconds"].empty() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test missing dependent option
 */
BOOST_AUTO_TEST_CASE( missing_dependent_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->depends("steps"), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::dependent_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK( vm["steps"].empty() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test reversed missing dependent option
 */
BOOST_AUTO_TEST_CASE( reversed_missing_dependent_option )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->depends("steps"), "")
        ("steps", po::value<uint64_t>(), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::dependent_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK( vm["steps"].empty() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test many dependent options
 */
BOOST_AUTO_TEST_CASE( many_dependent_options )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("seconds", po::value<double>(), "")
        ("time", po::value<double>()->depends("steps")->depends("seconds"), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296"
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::dependent_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["seconds"].empty() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test dependent options
 */
BOOST_AUTO_TEST_CASE( dependent_options )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("seconds", po::value<double>(), "")
        ("time", po::value<double>()->depends("steps")->depends("seconds"), "")
        ;
    array<char const*, 7> args = {{ "" //< argv[0]
      , "--steps", "4294967296"
      , "--seconds", "123.4567e-15"
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["seconds"].as<double>(), 123.4567e-15 );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test dependent options with defaulted values
 */
BOOST_AUTO_TEST_CASE( dependent_options_defaulted )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>()->default_value(4294967296LLU), "")
        ("time", po::value<double>()->depends("steps")->default_value(123.4567), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["steps"].defaulted() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( vm["time"].defaulted() );
}

/**
 * test dependent options with defaulted first value
 */
BOOST_AUTO_TEST_CASE( dependent_options_defaulted_first )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>()->default_value(4294967296LLU), "")
        ("time", po::value<double>()->depends("steps"), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::dependent_option>(parsed, vm, "time") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["steps"].defaulted() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( !vm["time"].defaulted() );
}

/**
 * test dependent options with defaulted second value
 */
BOOST_AUTO_TEST_CASE( dependent_options_defaulted_second )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->depends("steps")->default_value(123.4567), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--steps", "4294967296"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( !vm["steps"].defaulted() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( vm["time"].defaulted() );
}

/**
 * test non-existent dependent option
 */
BOOST_AUTO_TEST_CASE( non_existent_dependent_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("time", po::value<double>()->depends("step"), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296" //< MAX_UINT + 1LLU
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::unknown_option>(parsed, vm, "step") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test non-existent unused dependent option
 */
BOOST_AUTO_TEST_CASE( non_existent_unused_dependent_option )
{
    po::options_description desc;
    desc.add_options()
        ("steps", po::value<uint64_t>(), "")
        ("seconds", po::value<double>()->depends("step"), "")
        ("time", po::value<double>(), "")
        ;
    array<char const*, 5> args = {{ "" //< argv[0]
      , "--steps", "4294967296" //< MAX_UINT + 1LLU
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::unknown_option>(parsed, vm, "step") );
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["steps"].as<uint64_t>(), 4294967296LLU );
    BOOST_CHECK( vm["seconds"].empty() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test value store pointer
 */
BOOST_AUTO_TEST_CASE( value_store_pointer )
{
    po::options_description desc;
    double time = 0;
    desc.add_options()
        ("time", po::value<double>(&time), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( time, 123.4567 );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test default value
 */
BOOST_AUTO_TEST_CASE( default_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->default_value(123.4567), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( vm["time"].defaulted() );
}

/**
 * test overriden default value
 */
BOOST_AUTO_TEST_CASE( overriden_default_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->default_value(123.4567), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.45678"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_NE( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.45678 );
    BOOST_CHECK( !vm["time"].defaulted() );
}

/**
 * test implicit value
 */
BOOST_AUTO_TEST_CASE( implicit_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->implicit_value(123.4567), "")
        ;
    array<char const*, 2> args = {{ "" //< argv[0]
      , "--time"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK( !vm["time"].defaulted() );
}

/**
 * test empty implicit value
 */
BOOST_AUTO_TEST_CASE( empty_implicit_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->implicit_value(123.4567), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK( vm["time"].empty() );
}

/**
 * test overriden implicit value
 */
BOOST_AUTO_TEST_CASE( overriden_implicit_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->implicit_value(123.4567), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.45678"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_NE( vm["time"].as<double>(), 123.4567 );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.45678 );
    BOOST_CHECK( !vm["time"].defaulted() );
}

struct notify_exception
{
    int value;
    notify_exception(int value) : value(value) {}
};

void throw_notify(int value)
{
    throw notify_exception(value);
}

test_tools::predicate_result
notify_throws(po::variables_map& vm, int value)
{
    try {
        po::notify(vm);
    }
    catch (notify_exception const& e) {
        if (e.value != value) {
            test_tools::predicate_result res(false);
            res.message()
                << "Exception has wrong value [" << e.value
                << " != " << value << "]"
                ;
            return res;
        }
        return true;
    }
    test_tools::predicate_result res(false);
    res.message() << "No exception thrown";
    return res;
}

/**
 * test notifier
 */
BOOST_AUTO_TEST_CASE( notifier )
{
    po::options_description desc;
    desc.add_options()
        ("dimension", po::value<int>()->notifier(&throw_notify), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--dimension", "42"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK( notify_throws(vm, 42) );
    BOOST_CHECK_EQUAL( vm["dimension"].as<int>(), 42 );
}

#if BOOST_VERSION >= 104200

/**
 * test required value
 */
BOOST_AUTO_TEST_CASE( required_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->required(), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.45678"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.45678 );
}

/**
 * test missing required value
 */
BOOST_AUTO_TEST_CASE( missing_required_value )
{
    po::options_description desc;
    desc.add_options()
        ("time", po::value<double>()->required(), "")
        ;
    array<char const*, 1> args = {{ "" //< argv[0]
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    BOOST_CHECK_THROW( po::notify(vm), po::required_option );
    BOOST_CHECK( vm["time"].empty() );
}

#endif /* BOOST_VERSION >= 104200 */

/**
 * test multi-token value
 */
BOOST_AUTO_TEST_CASE( multitoken_value )
{
    po::options_description desc;
    desc.add_options()
        ("times", po::value<vector<double> >()->multitoken(), "")
        ("steps", po::value<vector<uint64_t> >()->multitoken(), "")
        ("time", po::value<double>(), "")
        ;
    array<char const*, 11> args = {{ "" //< argv[0]
      , "--times", "123.45678", "123.4567", "123.4567", "123.45679"
      , "--steps", "4294967296", "4294967296"
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >().size(), 4LU );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[0], 123.45678 );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[1], 123.4567 );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[2], 123.4567 );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[3], 123.45679 );
    BOOST_CHECK_EQUAL( vm["steps"].as<vector<uint64_t> >().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["steps"].as<vector<uint64_t> >()[0], 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["steps"].as<vector<uint64_t> >()[1], 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test composing value
 */
BOOST_AUTO_TEST_CASE( composing_value )
{
    po::options_description desc;
    desc.add_options()
        ("times", po::value<vector<double> >()->composing(), "")
        ("steps", po::value<vector<uint64_t> >()->composing(), "")
        ("time", po::value<double>(), "")
        ;
    array<char const*, 13> args = {{ "" //< argv[0]
      , "--times", "123.45678"
      , "--times", "123.4567"
      , "--times", "123.45679"
      , "--time", "123.4567"
      , "--steps", "4294967296"
      , "--steps", "4294967296"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >().size(), 3LU );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[0], 123.45678 );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[1], 123.4567 );
    BOOST_CHECK_EQUAL( vm["times"].as<vector<double> >()[2], 123.45679 );
    BOOST_CHECK_EQUAL( vm["steps"].as<vector<uint64_t> >().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["steps"].as<vector<uint64_t> >()[0], 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["steps"].as<vector<uint64_t> >()[1], 4294967296LLU );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

//
// zero_tokens is not tested, for details see:
//
// https://svn.boost.org/trac/boost/ticket/1132#comment:10
//
// > Therefore, 'zero_tokens' makes sense exclusively when the validator for
// > the type handles empty string or 'implied_value' is also specified.
//

/**
 * test true boolean value
 */
BOOST_AUTO_TEST_CASE( bool_switch_true )
{
    po::options_description desc;
    desc.add_options()
        ("debug", po::bool_switch(), "")
        ("time", po::value<double>(), "")
        ;
    array<char const*, 4> args = {{ "" //< argv[0]
      , "--debug"
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK( vm["debug"].as<bool>() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test false boolean value
 */
BOOST_AUTO_TEST_CASE( bool_switch_false )
{
    po::options_description desc;
    desc.add_options()
        ("debug", po::bool_switch(), "")
        ("time", po::value<double>(), "")
        ;
    array<char const*, 3> args = {{ "" //< argv[0]
      , "--time", "123.4567"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK( !vm["debug"].as<bool>() );
    BOOST_CHECK_EQUAL( vm["time"].as<double>(), 123.4567 );
}

/**
 * test untyped value
 */
BOOST_AUTO_TEST_CASE( untyped_value )
{
    po::options_description desc;
    desc.add_options()
        ("verbose", "")
        ;
    array<char const*, 2> args = {{ "" //< argv[0]
      , "--verbose"
    }};
    po::command_line_parser parser(args.size(), const_cast<char**>(&args.front()));
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("verbose"), 1LU );
}
