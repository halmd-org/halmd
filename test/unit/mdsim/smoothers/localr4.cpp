/*
 * Copyright © 2011 Felix Höfling
 * Copyright © 2012 Nicolas Höft
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE localr4
#include <boost/test/unit_test.hpp>

#ifndef HALMD_NO_CXX11

#include <array>
#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>
#include <limits>

#include <halmd/mdsim/host/potentials/lennard_jones.hpp>
#include <halmd/mdsim/smoothers/localr4.hpp>
#include <test/tools/ctest.hpp>

BOOST_AUTO_TEST_SUITE( host )

BOOST_AUTO_TEST_CASE( localr4 )
{
    typedef halmd::mdsim::host::potentials::lennard_jones<double> potential_type;
    typedef halmd::mdsim::smoothers::localr4<double> smooth_type;
    typedef potential_type::matrix_type matrix_type;

    float wca_cut = std::pow(2., 1 / 6.);
    // two particle types, i.e. binary mixture
    unsigned int ntype = 2;
    // cutoff at r_c=2.5σ
    matrix_type cutoff_array(ntype, ntype);
    cutoff_array <<=
        2.5, 2.
      , 2. , wca_cut;
    matrix_type epsilon_array(ntype, ntype);
    epsilon_array <<=
        1. , 0.5
      , 0.5, 0.25;
    matrix_type sigma_array(ntype, ntype);
    sigma_array <<=
        1., 2.
      , 2., 4.;
    // smoothing parameter
    double const h = 1. / 256;

    // construct potential module
    potential_type potential(ntype, ntype, cutoff_array, epsilon_array, sigma_array);
    smooth_type smooth(h);

    double const eps = std::numeric_limits<double>::epsilon();
    typedef std::array<double, 3> row_type;

    // expected results (r, fval, en_pot) for ε=1, σ=1, rc=2.5σ
    std::array<row_type, 5> const results_aa = {{
        {{0.5, 780287.9999895841, 16128.01631665645}}
      , {{2.4968, -0.01832978748170204, -0.00003892483911376598}}
      , {{2.491, -0.0176021892083538, -0.0003432476754602092}}
      , {{2.498, -0.00477727409418668, -5.029365367990809e-6}}
      , {{2.499, -0.0003331178245421451, -1.670172342312599e-7}}
    }};

    for (row_type const& row : results_aa) {
        double rr = std::pow(row[0], 2);
        double const rcut = potential.r_cut(0, 0);

        double fval, en_pot, hvir;
        // AA interaction
        boost::tie(fval, en_pot, hvir) = potential(rr, 0, 0);

        smooth(row[0], rcut, fval, en_pot);

        double const tolerance = 8 * eps * (1 + rcut / (rcut - row[0]));

        BOOST_CHECK_CLOSE_FRACTION(fval, row[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, row[2], tolerance);
    };

    // expected results (r, fval, en_pot) for ε=0.5, σ=2, rc=2σ
    std::array<row_type, 4> const results_ab = {{
        {{0.5, 1.61041612799762e9, 3.354624003070967e7}}
      , {{3.99, -0.01233402752055775, -0.0004475690049274307}}
      , {{3.997, -0.01167145018465577, -0.00003525113657243568}}
      , {{3.999, -0.0002422286662913668, -1.943663814679904e-7}}
    }};

    for (row_type const& row : results_ab) {
        double rr = std::pow(row[0], 2);
        double const rcut = potential.r_cut(0, 1);

        double fval, en_pot, hvir;
        // AB interaction
        boost::tie(fval, en_pot, hvir) = potential(rr, 0, 1);

        smooth(row[0], rcut, fval, en_pot);

        double const tolerance = 8 * eps * (1 + rcut / (rcut - row[0]));

        BOOST_CHECK_CLOSE_FRACTION(fval, row[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, row[2], tolerance);
    };

    // expected results (r, fval, en_pot) for ε=0.25, σ=4, rc=2^(1/6)σ
    std::array<row_type, 5> const results_bb = {{
        {{4., 0.3750000005473837, 0.2499999989890378}}
      , {{4.47, 0.004159744955159131, 0.0001811606617754403}}
      , {{4.48, 0.00205409992899728, 0.00004290963566710586}}
      , {{4.482, 0.001672906506430749, 0.00002622848763338342}}
      , {{4.487, 0.0003213380964310363, 8.01591174029315e-7}}
    }};

    for (row_type const& row : results_bb) {
        double rr = std::pow(row[0], 2);
        double fval, en_pot, hvir;
        // BB interaction
        boost::tie(fval, en_pot, hvir) = potential(rr, 1, 1);
        double const rcut = potential.r_cut(1, 1);

        smooth(row[0], rcut, fval, en_pot);
        double const tolerance = 9 * eps * (1 + rcut / (rcut - row[0]));

        BOOST_CHECK_CLOSE_FRACTION(fval, row[1], 2 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, row[2], tolerance);
    };
}

BOOST_AUTO_TEST_SUITE_END() // host

#endif /* ! HALMD_NO_CXX11 */
