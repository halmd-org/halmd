/*
 * Copyright © 2010-2012  Felix Höfling
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

#define BOOST_TEST_MODULE pow
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>
#include <cmath>
#include <limits>

#include <halmd/numeric/pow.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>

const double eps = std::numeric_limits<double>::epsilon();

using namespace halmd;

//
// test and benchmark fixed_pow() and halmd::pow() functions
//

BOOST_AUTO_TEST_CASE( test_fixed_pow )
{
    BOOST_CHECK_EQUAL(std::pow(2., 0), fixed_pow<0>(2.));
    BOOST_CHECK_EQUAL(std::pow(2., 1), fixed_pow<1>(2.));
    BOOST_CHECK_EQUAL(std::pow(2., 2), fixed_pow<2>(2.));
    BOOST_CHECK_EQUAL(std::pow(2., 6), fixed_pow<6>(2.));
    BOOST_CHECK_EQUAL(std::pow(2., 12), fixed_pow<12>(2.));
    BOOST_CHECK_EQUAL(std::pow(2., 15), fixed_pow<15>(2.));
    BOOST_CHECK_EQUAL(std::pow(2., 24), fixed_pow<24>(2.));

    BOOST_CHECK_CLOSE_FRACTION(fixed_pow<2>(std::sqrt(5.)), 5., eps);
    BOOST_CHECK_CLOSE_FRACTION(fixed_pow<12>(1.3), 23.298085122481, 11 * eps);  // 12-1 multiplications
}

BOOST_AUTO_TEST_CASE( test_halmd_pow )
{
    BOOST_CHECK_EQUAL(std::pow(2., 0), halmd::pow(2., 0));
    BOOST_CHECK_EQUAL(std::pow(2., 1), halmd::pow(2., 1));
    BOOST_CHECK_EQUAL(std::pow(2., 2), halmd::pow(2., 2));
    BOOST_CHECK_EQUAL(std::pow(2., 6), halmd::pow(2., 6));
    BOOST_CHECK_EQUAL(std::pow(2., 12), halmd::pow(2., 12));
    BOOST_CHECK_EQUAL(std::pow(2., 15), halmd::pow(2., 15));
    BOOST_CHECK_EQUAL(std::pow(2., 24), halmd::pow(2., 24));

    BOOST_CHECK_CLOSE_FRACTION(halmd::pow(std::sqrt(5.), 2), 5, eps);
    BOOST_CHECK_CLOSE_FRACTION(halmd::pow(1.3, 12), 23.298085122481, 11 * eps);
}

BOOST_AUTO_TEST_CASE( performance )
{
    boost::array<double, 5> elapsed;
    double a;
    unsigned n;

    // a const variable, even with known value, can not be used to
    // instantiate a template parameter
    unsigned const index = 12;
    enum { index_ = 12 };

    BOOST_TEST_MESSAGE("evaluation time of x^" << index << " in nanoseconds:");

    halmd::timer timer;
    a = 0;
    for (n=0; n < 100000000; n++) {
        a += fixed_pow<index_>((double)n);
    }
    elapsed[0] = timer.elapsed() / n;
    BOOST_TEST_MESSAGE("  fixed_pow: " << 1e9 * elapsed[0]);
    // the error bounds are a worst case scenario, the first term arises from the
    // error of exponentiation, and the second one from the summation of increasingly
    // large numbers; compare Sum[n^k, {n, 0, m}] and Sum[l^12, {n, 0, m}, {l, 0, n}].
    BOOST_CHECK_CLOSE_FRACTION(a, 7.692307192307702e102, (index_ + n / index_) * eps);

    timer.restart();
    a = 0;
    for (n=0; n < 100000000; n++) {
        a += halmd::pow((double)n, index);
    }
    elapsed[1] = timer.elapsed() / n;
    BOOST_TEST_MESSAGE("  halmd::pow: " << 1e9 * elapsed[1]);
    BOOST_CHECK_CLOSE_FRACTION(a, 7.692307192307702e102, (index + n / index) * eps);

    timer.restart();
    a = 0;
    for (n=0; n < 100000000; n++) {
        a += std::pow((double)n, index_);
    }
    elapsed[2] = timer.elapsed() / n;
    BOOST_TEST_MESSAGE("  std::pow (fixed): " << 1e9 * elapsed[2]);
    BOOST_CHECK_CLOSE_FRACTION(a, 7.692307192307702e102, (index_ + n / index_) * eps);

    timer.restart();
    a = 0;
    for (n=0; n < 1000000; n++) {
        a += std::pow((double)n, index); // seems to results in a call to libm
    }
    elapsed[3] = timer.elapsed() / n;
    BOOST_TEST_MESSAGE("  std::pow (variable): " << 1e9 * elapsed[3]);
    BOOST_CHECK_CLOSE_FRACTION(a, 7.692257692407692e76, (index + n / index) * eps);

    timer.restart();
    a = 0;
    for (n=0; n < 100000; n++) {
        a += std::pow((double)n, double(index));
    }
    elapsed[4] = timer.elapsed() / n;
    BOOST_TEST_MESSAGE("  std::pow (double): " << 1e9 * elapsed[4]);
    BOOST_CHECK_CLOSE_FRACTION(a, 7.691807702307692e63, (index + n / index) * eps);

    // With a recent GCC compiler (4.3, 4.4) and C++ standard library,
    // std::pow may in fact perform better than fixed_pow, so the original
    // elapsed time check has been downgraded to an informational message.

    BOOST_TEST_MESSAGE("  fixed_pow speed-up over std::pow with fixed exponent: " << elapsed[2] / elapsed[0]);
    BOOST_TEST_MESSAGE("  halmd::pow speed-up over std::pow with variable exponent: " << elapsed[3] / elapsed[1]);
}
