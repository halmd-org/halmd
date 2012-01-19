/*
 * Copyright © 2011  Felix Höfling
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

#define BOOST_TEST_MODULE pick_lattice_points
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/shared_ptr.hpp>
#include <iterator>
#include <numeric>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace halmd::algorithm::host;
using namespace std;

template <int dimension>
void pick_lattice_points()
{
    typedef fixed_vector<double, dimension> vector_type;
    typedef fixed_vector<unsigned int, dimension> index_type;

    double epsilon = 0.05;
    unsigned int max_count = 7;
    const vector_type unit_cell =
        (dimension == 3) ? list_of(2)(3)(5) : list_of(1)(1);

    vector<double> radii = list_of(.5)(1)(1.5)(1.91)(2)(2.7)(3)(3)(4)(5)(30);
    vector<unsigned int> count = (dimension == 3)
        ? list_of(0)(0)(0)(1)(1)(0)(2)(2)(1)(2)(max_count)
        : list_of(0)(2)(0)(2)(2)(1)(4)(4)(4)(6)(max_count);

    multimap<double, vector_type> lattice_points;

    // call with empty list of radii
    BOOST_MESSAGE("test empty range");
    vector<double> empty_range;
    pick_lattice_points_from_shell(
        empty_range.begin(), empty_range.end()
      , inserter(lattice_points, lattice_points.begin())
      , unit_cell, epsilon, max_count
    );
    BOOST_CHECK_EQUAL(lattice_points.size(), 0u);

    // call with zero box size
    BOOST_MESSAGE("test zero box size");
    pick_lattice_points_from_shell(
        radii.begin(), radii.end()
      , inserter(lattice_points, lattice_points.begin())
      , vector_type(0), epsilon, max_count
    );
    BOOST_CHECK_EQUAL(lattice_points.size(), 0u);

    // call with zero tolerance
    BOOST_MESSAGE("test zero tolerance");
    pick_lattice_points_from_shell(
        radii.begin(), radii.end()
      , inserter(lattice_points, lattice_points.begin())
      , unit_cell, 0., max_count //< tolerance must be a floating-point type!
    );
    BOOST_CHECK_EQUAL(lattice_points.size(), (dimension == 3) ? 11u : 18u);
    lattice_points.clear();

    // call with zero max_count
    BOOST_MESSAGE("test zero max_count");
    pick_lattice_points_from_shell(
        radii.begin(), radii.end()
      , inserter(lattice_points, lattice_points.begin())
      , unit_cell, epsilon, 0
    );
    BOOST_CHECK_EQUAL(lattice_points.size(), 0u);

    // construct lattice points
    BOOST_MESSAGE("construct lattice points");
    pick_lattice_points_from_shell(
        radii.begin(), radii.end()
      , inserter(lattice_points, lattice_points.begin())
      , unit_cell, epsilon, max_count
    );
    BOOST_CHECK_EQUAL(lattice_points.size(), (dimension == 3) ? 14u : 28u);  // not equal to sum(count)

    // check conditions and counts on constructed lattice points
    for (unsigned int i = 0; i < radii.size(); ++i) {
        double r = radii[i];
        // check total count per shell
        unsigned int c = lattice_points.count(r);
        BOOST_CHECK_MESSAGE(c == count[i], c << " != " << count[i] << " for r = " << r);
        BOOST_CHECK_MESSAGE(c <= max_count, "count(r = " << r << ") = " << c);

        typedef typename multimap<double, vector_type>::const_iterator iterator_type;
        typedef pair<iterator_type, iterator_type> range_type;
        unsigned int sum = 0;
        for (range_type shell = lattice_points.equal_range(r); shell.first != shell.second; ++shell.first) {
            // check that distance to origin is within the tolerance
            vector_type const& point = shell.first->second;
            BOOST_CHECK_SMALL(norm_2(point) / r - 1, epsilon);

            // check ascending sum of Miller indices
            index_type hkl = static_cast<index_type>(round(element_div(point, unit_cell)));
            hkl /= greatest_common_divisor(hkl);
            unsigned int sum_ = accumulate(hkl.begin(), hkl.end(), 0u, plus<unsigned int>());
            BOOST_CHECK_MESSAGE(sum_ >= sum, "incorrect order of lattice points");
            sum_ = sum;
        }
    }
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test::framework;

    master_test_suite().add(BOOST_TEST_CASE( &pick_lattice_points<2> ));
    master_test_suite().add(BOOST_TEST_CASE( &pick_lattice_points<3> ));
}
