/*
 * Copyright © 2012      Nicolas Höft
 * Copyright © 2011-2013 Felix Höfling
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

#define BOOST_TEST_MODULE local_r4
#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <cmath>
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/mdsim/host/potentials/lennard_jones.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/lennard_jones.hpp>
# include <test/unit/mdsim/potentials/gpu/neighbour_chain.hpp>
# include <test/tools/cuda.hpp>
#endif
#include <test/tools/ctest.hpp>

BOOST_AUTO_TEST_SUITE( host )

BOOST_AUTO_TEST_CASE( local_r4 )
{
    typedef halmd::mdsim::host::potentials::lennard_jones<double> potential_type;
    typedef halmd::mdsim::forces::trunc::local_r4<double> trunc_type;
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
    potential_type potential(cutoff_array, epsilon_array, sigma_array);
    trunc_type trunc(h);

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

        double fval, en_pot;
        // AA interaction
        boost::tie(fval, en_pot) = potential(rr, 0, 0);

        trunc(row[0], rcut, fval, en_pot);

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

        double fval, en_pot;
        // AB interaction
        boost::tie(fval, en_pot) = potential(rr, 0, 1);

        trunc(row[0], rcut, fval, en_pot);

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
        double fval, en_pot;
        // BB interaction
        boost::tie(fval, en_pot) = potential(rr, 1, 1);
        double const rcut = potential.r_cut(1, 1);

        trunc(row[0], rcut, fval, en_pot);
        double const tolerance = 9 * eps * (1 + rcut / (rcut - row[0]));

        BOOST_CHECK_CLOSE_FRACTION(fval, row[1], 2 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, row[2], tolerance);
    };
}

BOOST_AUTO_TEST_SUITE_END() // host

#ifdef HALMD_WITH_GPU

BOOST_AUTO_TEST_SUITE( gpu )

template <typename float_type>
struct test_local_r4
{
    enum { dimension = 2 };

    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::forces::trunc::local_r4<float_type> trunc_type;
    typedef halmd::mdsim::forces::trunc::local_r4<double> host_trunc_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::potentials::lennard_jones<float_type> potential_type;
    typedef halmd::mdsim::host::potentials::lennard_jones<double> host_potential_type;
    typedef halmd::mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type, trunc_type> force_type;
    typedef neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    std::vector<unsigned int> npart_list;
    std::shared_ptr<box_type> box;
    std::shared_ptr<potential_type> potential;
    std::shared_ptr<trunc_type> trunc;
    std::shared_ptr<host_trunc_type> host_trunc;
    std::shared_ptr<force_type> force;
    std::shared_ptr<neighbour_type> neighbour;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<host_potential_type> host_potential;

    test_local_r4();
    void test();
};

template <typename float_type>
void test_local_r4<float_type>::test()
{
    // place particles along the x-axis within one half of the box,
    // put every second particle at the origin
    unsigned int npart = particle->nparticle();
    vector_type dx(0);
    dx[0] = box->edges()(0, 0) / npart / 2;

    std::vector<vector_type> r_list(particle->nparticle());
    std::vector<unsigned int> species(particle->nparticle());
    for (unsigned int k = 0; k < r_list.size(); ++k) {
        r_list[k] = (k % 2) ? k * dx : vector_type(0);
        species[k] = (k < npart_list[0]) ? 0U : 1U;  // set particle type for a binary mixture
    }
    BOOST_CHECK( set_position(*particle, r_list.begin()) == r_list.end() );
    BOOST_CHECK( set_species(*particle, species.begin()) == species.end() );

    // read forces and other stuff from device
    std::vector<float> en_pot(particle->nparticle());
    BOOST_CHECK( get_potential_energy(*particle, en_pot.begin()) == en_pot.end() );

    std::vector<vector_type> f_list(particle->nparticle());
    BOOST_CHECK( get_force(*particle, f_list.begin()) == f_list.end() );

    float_type const eps = std::numeric_limits<float_type>::epsilon();

    for (unsigned int i = 0; i < npart; ++i) {
        unsigned int type1 = species[i];
        unsigned int type2 = species[(i + 1) % npart];
        vector_type r = r_list[i] - r_list[(i + 1) % npart];
        vector_type f = f_list[i];

        // reference values from host module
        double fval, en_pot_;
        double rr = inner_prod(r, r);
        boost::tie(fval, en_pot_) = (*host_potential)(rr, type1, type2);

        if (rr < host_potential->rr_cut(type1, type2)) {
            double rcut = host_potential->r_cut(type1, type2);
            (*host_trunc)(std::sqrt(rr), rcut, fval, en_pot_);
            // the GPU force module stores only a fraction of these values
            en_pot_ /= 2;

            // the first term is from the smoothing, the second from the potential
            // (see lennard_jones.cpp from unit tests)
            float_type const tolerance = 8 * eps * (1 + rcut/(rcut - std::sqrt(rr))) + 10 * eps;

            BOOST_CHECK_SMALL(norm_inf(fval * r - f), std::max(norm_inf(fval * r), float_type(1)) * tolerance * 2);
            BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 2 * tolerance);
        }
        else {
            // when the distance is greater than the cutoff
            // set the force and the pair potential to zero
            fval = en_pot_ = 0;
            BOOST_CHECK_SMALL(norm_inf(f), eps);
            BOOST_CHECK_SMALL(en_pot[i], eps);
        }
    }
}

template <typename float_type>
test_local_r4<float_type>::test_local_r4()
{
    typedef typename potential_type::matrix_type matrix_type;

    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart_list = {500, 300};
    // the box length should not be greater than 2*r_c
    float box_length = 10;
    unsigned int const dimension = box_type::vector_type::static_size;
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }
    // smoothing parameter
    double const h = 1. / 256;

    float_type wca_cut = std::pow(2., 1. / 6);
    matrix_type cutoff_array(2, 2);
    cutoff_array <<=
        2.5, 2.
      , 2. , wca_cut;
    matrix_type epsilon_array(2, 2);
    epsilon_array <<=
        1., 2.
      , 2., 0.25;
    matrix_type sigma_array(2, 2);
    sigma_array <<=
        1., 2.
      , 2., 4.;

    // create modules
    particle = std::make_shared<particle_type>(accumulate(npart_list.begin(), npart_list.end(), 0), npart_list.size());
    box = std::make_shared<box_type>(edges);
    potential = std::make_shared<potential_type>(cutoff_array, epsilon_array, sigma_array);
    host_potential = std::make_shared<host_potential_type>(cutoff_array, epsilon_array, sigma_array);
    neighbour = std::make_shared<neighbour_type>(particle);
    trunc = std::make_shared<trunc_type>(h);
    host_trunc = std::make_shared<host_trunc_type>(h);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour, trunc);
    particle->on_prepend_force([=](){force->check_cache();});
    particle->on_force([=](){force->apply();});
}

BOOST_FIXTURE_TEST_CASE( local_r4, set_cuda_device ) {
    test_local_r4<float>().test();
}

BOOST_AUTO_TEST_SUITE_END() // gpu

#endif /** HALMD_WITH_GPU */
