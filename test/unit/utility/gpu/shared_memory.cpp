/*
 * Copyright © 2012  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE shared_memory
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>

#include <halmd/utility/gpu/shared_memory.hpp>
#include <test/tools/ctest.hpp>

using namespace halmd;

BOOST_AUTO_TEST_CASE( test_power_of_two_max_threads )
{
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<32>::value), 32 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<33>::value), 32 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<42>::value), 32 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<63>::value), 32 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<64>::value), 64 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<65>::value), 64 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<127>::value), 64 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<128>::value), 128 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<129>::value), 128 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<255>::value), 128 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<256>::value), 256 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<257>::value), 256 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<511>::value), 256 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<512>::value), 512 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<513>::value), 512 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<1023>::value), 512 );
#if HALMD_GPU_ARCH >= 200
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<1024>::value), 1024 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<1025>::value), 1024 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<2047>::value), 1024 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<2048>::value), 1024 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<2049>::value), 1024 );
#else
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<1024>::value), 512 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<1025>::value), 512 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<2047>::value), 512 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<2048>::value), 512 );
    BOOST_CHECK_EQUAL( int(power_of_two_max_threads<2049>::value), 512 );
#endif
}

typedef boost::array<float, 7> float7;
typedef boost::array<double, 4> double4;
typedef boost::array<double, 8> double8;

BOOST_AUTO_TEST_CASE( test_shared_memory_max_threads )
{
#if HALMD_GPU_ARCH >= 200
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<char>::value), 1024 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<short>::value), 1024 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<float>::value), 1024 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<double>::value), 1024 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<float7>::value), 1024 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<double4>::value), 1024 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<double8>::value), 512 );
#else
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<char>::value), 512 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<short>::value), 512 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<float>::value), 512 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<double>::value), 512 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<float7>::value), 512 );
    // with PTX 1.3 and older, accounts for shared memory used by kernel arguments
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<double4>::value), 256 );
    BOOST_CHECK_EQUAL( int(shared_memory_max_threads<double8>::value), 128 );
#endif
}
