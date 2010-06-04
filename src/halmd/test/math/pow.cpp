/* test and benchmark math::pow() function
 *
 * Copyright (C) 2010  Felix Höfling
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_math_pow
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <limits>

#include <halmd/math/pow.hpp>
#include <halmd/util/timer.hpp>

const double eps = std::numeric_limits<double>::epsilon();

using namespace halmd;

BOOST_AUTO_TEST_CASE( correctness )
{
    BOOST_CHECK_EQUAL(std::pow(2, 0), math::pow<0>(2));
    BOOST_CHECK_EQUAL(std::pow(2, 1), math::pow<1>(2));
    BOOST_CHECK_EQUAL(std::pow(2, 2), math::pow<2>(2));
    BOOST_CHECK_EQUAL(std::pow(2, 6), math::pow<6>(2));
    BOOST_CHECK_EQUAL(std::pow(2, 12), math::pow<12>(2));
    BOOST_CHECK_EQUAL(std::pow(2, 15), math::pow<15>(2));
    BOOST_CHECK_EQUAL(std::pow(2, 24), math::pow<24>(2));

    BOOST_CHECK_CLOSE_FRACTION(math::pow<2>(std::sqrt(5)), 5, eps);
    BOOST_CHECK_CLOSE_FRACTION(math::pow<12>(1.3), std::pow(1.3, 12), eps);
}

BOOST_AUTO_TEST_CASE( performance )
{
    high_resolution_timer timer[6];
    unsigned n;
    #define index 12

    BOOST_TEST_MESSAGE("evaluation time of x^" << index << " in nanoseconds:");
    timer[0].record();
    double a = 0;
    for (n=0; n < 10000000; n++) {
        a += math::pow<index>((double)n);
    }
    timer[1].record();
    BOOST_TEST_MESSAGE("  math::pow: " << 1e9 * (timer[1] - timer[0]) / n);

    timer[2].record();
    double b = 0;
    for (n=0; n < 10000000; n++) {
        b += std::pow((double)n, index);
    }
    timer[3].record();
    BOOST_TEST_MESSAGE("  std::pow: " << 1e9 * (timer[3] - timer[2]) / n);

    timer[4].record();
    double c = 0;
    for (n=0; n < 100000; n++) {
        c += std::pow((double)n, double(index));
    }
    timer[5].record();
    BOOST_TEST_MESSAGE("  std::pow (double): " << 1e9 * (timer[5] - timer[4]) / n);

    BOOST_CHECK_CLOSE_FRACTION(a, b, 2 * eps);
    BOOST_CHECK_MESSAGE(timer[1] - timer[0] < .9 * (timer[3] - timer[2]),
                        "improvement of math::pow is less than 20% in speed");
}

