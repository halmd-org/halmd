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

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/utility/wavevectors.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace std;

using observables::utility::construct_wavevector_shells;

BOOST_AUTO_TEST_CASE( wavevectors )
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

    typedef fixed_vector<double, 3> vector_type;
    typedef fixed_vector<unsigned int, 3> index_type;

    vector<double> wavenumbers = list_of(0.3)(0.7)(1.0)(1.5)(2.0)(25.0);
    const vector_type box_length = list_of(10.)(10.)(20.);
    double epsilon = 0.1;
    unsigned int miller_max = 4;

    multimap<double, vector_type> wavevectors =
        construct_wavevector_shells(wavenumbers, box_length, epsilon, miller_max);

    // check conditions on constructed wavevectors
    BOOST_FOREACH (double q, wavenumbers) {
        typedef multimap<double, vector_type>::const_iterator iterator_type;
        typedef pair<iterator_type, iterator_type> range_type;
        for (range_type shell = wavevectors.equal_range(q); shell.first != shell.second; ++shell.first) {
            vector_type const& q_vector = shell.first->second;
            const vector_type q_basis = element_div(vector_type(2 * M_PI), box_length);
            index_type hkl = static_cast<index_type>(round(element_div(q_vector, q_basis)));
            hkl /= greatest_common_divisor(hkl);
            BOOST_CHECK_SMALL(norm_2(q_vector) / q - 1, epsilon);
            BOOST_CHECK(*max_element(hkl.begin(), hkl.end()) <= miller_max);
        }
    }
}
