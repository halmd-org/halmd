/*
 * Copyright © 2010  Felix Höfling
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

#define BOOST_TEST_MODULE pow
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>
#include <cmath>
#include <limits>

#include <halmd/numeric/pow.hpp>
#include <halmd/utility/timer.hpp>

const double eps = std::numeric_limits<double>::epsilon();

using namespace halmd;

//
// test and benchmark fixed_pow() function
//

BOOST_AUTO_TEST_CASE( correctness )
{
    BOOST_CHECK_EQUAL(std::pow(2., 0), fixed_pow<0>(2));
    BOOST_CHECK_EQUAL(std::pow(2., 1), fixed_pow<1>(2));
    BOOST_CHECK_EQUAL(std::pow(2., 2), fixed_pow<2>(2));
    BOOST_CHECK_EQUAL(std::pow(2., 6), fixed_pow<6>(2));
    BOOST_CHECK_EQUAL(std::pow(2., 12), fixed_pow<12>(2));
    BOOST_CHECK_EQUAL(std::pow(2., 15), fixed_pow<15>(2));
    BOOST_CHECK_EQUAL(std::pow(2., 24), fixed_pow<24>(2));

    BOOST_CHECK_CLOSE_FRACTION(fixed_pow<2>(std::sqrt(5.)), 5, eps);
    BOOST_CHECK_CLOSE_FRACTION(fixed_pow<12>(1.3), std::pow(1.3, 12), eps);
}

BOOST_AUTO_TEST_CASE( performance )
{
    boost::array<double, 3> elapsed;
    unsigned n;
    #define index 12

    BOOST_TEST_MESSAGE("evaluation time of x^" << index << " in nanoseconds:");
    halmd::timer timer;
    double a = 0;
    for (n=0; n < 10000000; n++) {
        a += fixed_pow<index>((double)n);
    }
    elapsed[0] = timer.elapsed();
    BOOST_TEST_MESSAGE("  fixed_pow: " << 1e9 * elapsed[0] / n);

    timer.restart();
    double b = 0;
    for (n=0; n < 10000000; n++) {
        b += std::pow((double)n, index);
    }
    elapsed[1] = timer.elapsed();
    BOOST_TEST_MESSAGE("  std::pow: " << 1e9 * elapsed[1] / n);

    timer.restart();
    double c = 0;
    for (n=0; n < 100000; n++) {
        c += std::pow((double)n, double(index));
    }
    elapsed[2] = timer.elapsed();
    BOOST_TEST_MESSAGE("  std::pow (double): " << 1e9 * elapsed[2] / n);

    // use the results in some way in order to avoid an overoptimization of the loops
    BOOST_TEST_MESSAGE("  results: " << a << " " << b << " " << c);

    BOOST_CHECK_CLOSE_FRACTION(a, b, 2 * eps);

    // With a recent GCC compiler (4.3, 4.4) and C++ standard library,
    // std::pow may in fact perform better than fixed_pow, so the original
    // elapsed time check has been downgraded to an informational message.

    BOOST_TEST_MESSAGE("  fixed_pow speed-up over std::pow: "
                       << 100 * (1 - elapsed[0] / elapsed[1]) << "%");
}
