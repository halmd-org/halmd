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

#define BOOST_TEST_MODULE test_accumulator
#include <boost/test/unit_test.hpp>

#include <halmd/numeric/accumulator.hpp>

BOOST_AUTO_TEST_CASE( test_accumulator )
{
    halmd::accumulator<double> a;

    for (unsigned i=0; i <= 1; i++) {
        for (unsigned j=0; j < 10; j++) {
            a(j);
        }

        BOOST_CHECK_EQUAL(count(a), 10u);
        BOOST_CHECK_CLOSE_FRACTION(mean(a), 4.5, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(variance(a), 8.25, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(sigma(a), 2.8722813232690143, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(error_of_mean(a), 0.9574271077563381, 1e-14);

        a = halmd::accumulator<double>();  // check value initialization
    }
}
