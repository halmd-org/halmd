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

#define BOOST_TEST_MODULE algorithm_host_pick_lattice_points
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <iterator>
#include <numeric>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace halmd::algorithm::host;
using namespace std;

template <int dimension>
void pick_lattice_points()
{
    // enable logging to console
    shared_ptr<logger> log(new logger);
    log->log_to_console(
#ifdef NDEBUG
        logger::warning
#else
        logger::trace
#endif
    );

    typedef fixed_vector<double, dimension> vector_type;
    typedef fixed_vector<unsigned int, dimension> index_type;

    vector<double> radii = list_of(0.3)(0.7)(1.0)(1.5)(2.0)(25.0);
    const vector_type unit_cell =
        (dimension == 3) ? list_of(.5)(.5)(.2) : list_of(.3)(.3);
    double epsilon = 0.05;
    unsigned int max_count = 10;

    // construct lattice points
    multimap<double, vector_type> lattice_points;
    pick_lattice_points_from_shell(
        radii.begin(), radii.end()
      , inserter(lattice_points, lattice_points.begin())
      , unit_cell, epsilon, max_count
    );

    // check conditions on constructed lattice points
    BOOST_FOREACH (double r, radii) {
        // check total count per shell
        unsigned count = lattice_points.count(r);
        BOOST_CHECK(count <= max_count);
        if (!count) {
            BOOST_TEST_MESSAGE("No lattice points compatible with r ≈ " << r);
        }

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

static void __attribute__((constructor)) init_unit_test_suite()
{
    using namespace boost::unit_test::framework;

    master_test_suite().add(BOOST_TEST_CASE( &pick_lattice_points<2> ));
    master_test_suite().add(BOOST_TEST_CASE( &pick_lattice_points<3> ));
}
