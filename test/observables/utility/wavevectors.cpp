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

#define BOOST_TEST_MODULE observables_utility_wavevectors
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <numeric>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/utility/wavevectors.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace std;

using observables::utility::construct_wavevector_shells;

template <int dimension>
void wavevectors()
{
    // enable logging to console
    shared_ptr<logger> log(new logger);
    log->log_to_console(
#ifdef NDEBUG
        logger::warning
#else
        logger::debug
#endif
    );

    typedef fixed_vector<double, dimension> vector_type;
    typedef fixed_vector<unsigned int, dimension> index_type;

    vector<double> wavenumbers = list_of(0.3)(0.7)(1.0)(1.5)(2.0)(25.0);
    const vector_type box_length =
        (dimension == 3) ? list_of(10.)(10.)(20.) : list_of(20.)(20.);
    double epsilon = 0.05;
    unsigned int max_count = 10;

    multimap<double, vector_type> wavevectors =
        construct_wavevector_shells(wavenumbers, box_length, epsilon, max_count);

    // check conditions on constructed wavevectors
    BOOST_FOREACH (double q, wavenumbers) {
        // check total count per shell
        BOOST_CHECK(wavevectors.count(q) <= max_count);
        typedef typename multimap<double, vector_type>::const_iterator iterator_type;
        typedef pair<iterator_type, iterator_type> range_type;
        unsigned int sum = 0;
        for (range_type shell = wavevectors.equal_range(q); shell.first != shell.second; ++shell.first) {
            // check that magnitude is within the tolerance
            vector_type const& q_vector = shell.first->second;
            BOOST_CHECK_SMALL(norm_2(q_vector) / q - 1, epsilon);

            // check ascending sum of Miller indices
            const vector_type q_basis = element_div(vector_type(2 * M_PI), box_length);
            index_type hkl = static_cast<index_type>(round(element_div(q_vector, q_basis)));
            hkl /= greatest_common_divisor(hkl);
            unsigned int sum_ = accumulate(hkl.begin(), hkl.end(), 0u, plus<unsigned int>());
            BOOST_CHECK_MESSAGE(sum_ >= sum, "incorrect order of wavevectors");
            sum_ = sum;
        }
    }
}

static void __attribute__((constructor)) init_unit_test_suite()
{
    using namespace boost::unit_test::framework;

    master_test_suite().add(BOOST_TEST_CASE( &wavevectors<2> ));
    master_test_suite().add(BOOST_TEST_CASE( &wavevectors<3> ));
}
