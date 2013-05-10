/*
 * Copyright © 2010  Felix Höfling
 * Copyright © 2013  Nicolas Höft
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

#define BOOST_TEST_MODULE accumulator
#include <boost/test/unit_test.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/tools/ctest.hpp>

BOOST_AUTO_TEST_CASE( test_accumulator )
{
    halmd::accumulator<double> a;

    for (unsigned int i = 0; i <= 1; ++i) {
        for (unsigned int j = 0; j < 10; ++j) {
            a(j);
        }

        BOOST_CHECK_EQUAL(count(a), 10u);
        BOOST_CHECK_CLOSE_FRACTION(mean(a), 4.5, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(variance(a), 8.25, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(sigma(a), 2.8722813232690143, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(error_of_mean(a), 0.9574271077563381, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(sum(a), 45., 1e-14);

        a = halmd::accumulator<double>();  // check value initialization
    }

    unsigned int constexpr dim = 2;
    halmd::accumulator<halmd::fixed_vector<double, dim>> av;
    halmd::accumulator<halmd::fixed_vector<double, dim>> av2;

    for (unsigned i = 0; i <= 1; ++i) {
        for (unsigned j = 0; j < 10; ++j) {
            halmd::fixed_vector<double, dim> v;
            for (unsigned int d = 0; d < dim; ++d) {
                v[d] = j;
            }
            av(v);
            av2(v);
        }

        BOOST_CHECK_EQUAL(count(av), 10u);
        for (unsigned int d = 0; d < dim; ++d) {
            BOOST_CHECK_CLOSE_FRACTION(mean(av)[d], 4.5, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(variance(av)[d], 8.25, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(sigma(av)[d], 2.8722813232690143, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(error_of_mean(av)[d], 0.9574271077563381, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(sum(av)[d], 45., 1e-14);
        }
        // test functor for accumulators
        av2(av);
        BOOST_CHECK_EQUAL(count(av2), 20u);
        for (unsigned int d = 0; d < dim; ++d) {
            BOOST_CHECK_CLOSE_FRACTION(mean(av2)[d], 4.5, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(variance(av2)[d], 8.25, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(sigma(av2)[d], 2.8722813232690143, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(error_of_mean(av2)[d], 0.65894652766046917965, 1e-14);
            BOOST_CHECK_CLOSE_FRACTION(sum(av2)[d], 90., 1e-14);
        }

        av = halmd::accumulator<halmd::fixed_vector<double, dim>>();  // check value initialization
        av2 = halmd::accumulator<halmd::fixed_vector<double, dim>>();
    }

}
