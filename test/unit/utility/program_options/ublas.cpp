/*
 * Copyright © 2011  Peter Colberg
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

#define BOOST_TEST_MODULE program_options_ublas
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <limits>

#include <halmd/utility/program_options/ublas.hpp>
#include <test/unit/utility/program_options/predicates.hpp>

using namespace boost;
using namespace boost::assign;
using namespace std;

namespace po = boost::program_options;
namespace ublas = boost::numeric::ublas;

BOOST_AUTO_TEST_CASE( matrix_1_by_1 )
{
    float const eps = numeric_limits<float>::epsilon();
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::matrix<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("3.141592653589793");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("epsilon"), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size1(), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size2(), 1LU );
    BOOST_CHECK_CLOSE_FRACTION( vm["epsilon"].as<ublas::matrix<float> >()(0, 0), M_PI, eps );
    BOOST_MESSAGE( vm["epsilon"].as<ublas::matrix<float> >() );
}

BOOST_AUTO_TEST_CASE( matrix_1_by_2 )
{
    float const eps = numeric_limits<float>::epsilon();
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::matrix<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("3.141592653589793,2.718281828459045");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("epsilon"), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size1(), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size2(), 2LU );
    BOOST_CHECK_CLOSE_FRACTION( vm["epsilon"].as<ublas::matrix<float> >()(0, 0), M_PI, eps );
    BOOST_CHECK_CLOSE_FRACTION( vm["epsilon"].as<ublas::matrix<float> >()(0, 1), M_E, eps );
    BOOST_MESSAGE( vm["epsilon"].as<ublas::matrix<float> >() );
}

BOOST_AUTO_TEST_CASE( matrix_2_by_1 )
{
    float const eps = numeric_limits<float>::epsilon();
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::matrix<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("3.141592653589793:2.718281828459045");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("epsilon"), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size1(), 2LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size2(), 1LU );
    BOOST_CHECK_CLOSE_FRACTION( vm["epsilon"].as<ublas::matrix<float> >()(0, 0), M_PI, eps );
    BOOST_CHECK_CLOSE_FRACTION( vm["epsilon"].as<ublas::matrix<float> >()(1, 0), M_E, eps );
    BOOST_MESSAGE( vm["epsilon"].as<ublas::matrix<float> >() );
}

BOOST_AUTO_TEST_CASE( matrix_2_by_2 )
{
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::matrix<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("1.1,1.2:2.1,2.2");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("epsilon"), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size1(), 2LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >().size2(), 2LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >()(0, 0), 1.1f );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >()(0, 1), 1.2f );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >()(1, 0), 2.1f );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::matrix<float> >()(1, 1), 2.2f );
    BOOST_MESSAGE( vm["epsilon"].as<ublas::matrix<float> >() );
}

BOOST_AUTO_TEST_CASE( matrix_throws_invalid_options_value )
{
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::matrix<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("1.1,1.2:2.1,2.2:2.2");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::invalid_option_value>(parsed, vm, "epsilon") );
}

BOOST_AUTO_TEST_CASE( vector_1 )
{
    float const eps = numeric_limits<float>::epsilon();
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::vector<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("3.141592653589793");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("epsilon"), 1LU );
    BOOST_CHECK_EQUAL( vm["epsilon"].as<ublas::vector<float> >().size(), 1LU );
    BOOST_CHECK_CLOSE_FRACTION( vm["epsilon"].as<ublas::vector<float> >()(0), M_PI, eps );
    BOOST_MESSAGE( vm["epsilon"].as<ublas::vector<float> >() );
}

BOOST_AUTO_TEST_CASE( vector_2 )
{
    po::options_description desc;
    desc.add_options()
        ("diameter", po::value<ublas::vector<float> >(), "")
        ;
    vector<string> args = list_of("--diameter")("1.1,1.2");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);
    BOOST_CHECK_EQUAL( vm.count("diameter"), 1LU );
    BOOST_CHECK_EQUAL( vm["diameter"].as<ublas::vector<float> >().size(), 2LU );
    BOOST_CHECK_EQUAL( vm["diameter"].as<ublas::vector<float> >()(0), 1.1f );
    BOOST_CHECK_EQUAL( vm["diameter"].as<ublas::vector<float> >()(1), 1.2f );
    BOOST_MESSAGE( vm["diameter"].as<ublas::vector<float> >() );
}

BOOST_AUTO_TEST_CASE( vector_throws_invalid_options_value )
{
    po::options_description desc;
    desc.add_options()
        ("epsilon", po::value<ublas::vector<float> >(), "")
        ;
    vector<string> args = list_of("--epsilon")("1.1:1.2");
    po::command_line_parser parser(args);
    po::parsed_options parsed(parser.options(desc).run());
    po::variables_map vm;
    BOOST_CHECK( store_throws_option<po::invalid_option_value>(parsed, vm, "epsilon") );
}
