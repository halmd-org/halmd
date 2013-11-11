/*
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

#define BOOST_TEST_MODULE lennard_jones_linear
#include <boost/test/unit_test.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/lennard_jones_linear.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/lennard_jones_linear.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/gpu/neighbour_chain.hpp>
#endif
#include <test/tools/ctest.hpp>

#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <numeric>

/** test Lennard-Jones potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

BOOST_AUTO_TEST_CASE( lennard_jones_linear_host )
{
    typedef halmd::mdsim::host::potentials::lennard_jones_linear<double> potential_type;
    typedef potential_type::matrix_type matrix_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    matrix_type cutoff_array(ntype, ntype);
    matrix_type epsilon_array(ntype, ntype);
    matrix_type sigma_array(ntype, ntype);
    cutoff_array <<= 5., 5., 5., 5.;
    epsilon_array <<= 1., .5, .5, .25;
    sigma_array <<= 1., 2., 2., 4.;

    // construct module
    potential_type potential(cutoff_array, epsilon_array, sigma_array);

    // test paramters
    matrix_type epsilon = potential.epsilon();
    BOOST_CHECK(epsilon(0, 0) == epsilon_array(0, 0));
    BOOST_CHECK(epsilon(0, 1) == epsilon_array(0, 1));
    BOOST_CHECK(epsilon(1, 0) == epsilon_array(1, 0));
    BOOST_CHECK(epsilon(1, 1) == epsilon_array(1, 1));

    matrix_type sigma = potential.sigma();
    BOOST_CHECK(sigma(0, 0) == sigma_array(0, 0));
    BOOST_CHECK(sigma(0, 1) == sigma_array(0, 1));
    BOOST_CHECK(sigma(1, 0) == sigma_array(1, 0));
    BOOST_CHECK(sigma(1, 1) == sigma_array(1, 1));

    // evaluate some points of potential and force
    typedef std::array<double, 3> array_type;
    const double tolerance = 5 * std::numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot) for ε=1, σ=1, rc=5σ
    std::array<array_type, 5> results_aa = {{
        {{0.2, 2.929593750000015e11, 9.765000000017304e8}}
      , {{0.5, 780288.0006143214, 16128.00163820667}}
      , {{1., 24.0003071606784, 0.0014846263296}}
      , {{2., -0.0906667321608, -0.0603459718488}}
      , {{10., 0.00003047606832, -0.001283819772}}
    }};

    for (array_type const& a : results_aa) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        boost::tie(fval, en_pot) = potential(rr, 0, 0);  // interaction AA
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction AB: ε=.5, σ=2, rc=5σ
    std::array<array_type, 5> results_ab = {{
        {{0.2, 5.999997e14, 1.999998000000001e12}}
      , {{0.5, 1.610416128000154e9, 3.35462400008575e7}}
      , {{1., 97536.00007679017, 8064.000819103334}}
      , {{2., 3.0000383950848, 0.0007423131648}}
      , {{10., 0, 0}}
    }};

    for (array_type const& a : results_ab) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        boost::tie(fval, en_pot) = potential(rr, 0, 1);  // interaction AB
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction BB: ε=.25, σ=4, rc=5σ
    std::array<array_type, 5> results_bb = {{
        {{0.2, 1.2287999904e18, 4.095999936e15}}
      , {{0.5, 3.298528591872e12, 6.871921459200044e10}}
      , {{1., 2.013020160000192e8, 1.677312000042875e7}}
      , {{2., 12192.00000959877, 4032.000409551667}}
      , {{10., -0.00024182697984, -0.003823251456}}
    }};

    for (array_type const& a : results_bb) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        boost::tie(fval, en_pot) = potential(rr, 1, 1);  // interaction BB
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };
}

#ifdef HALMD_WITH_GPU

template <typename float_type>
struct lennard_jones_linear
{
    enum { dimension = 2 };

    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::potentials::lennard_jones_linear<float_type> potential_type;
    typedef halmd::mdsim::host::potentials::lennard_jones_linear<double> host_potential_type;
    typedef halmd::mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    std::shared_ptr<box_type> box;
    std::shared_ptr<potential_type> potential;
    std::shared_ptr<force_type> force;
    std::shared_ptr<neighbour_type> neighbour;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<host_potential_type> host_potential;
    std::vector<unsigned int> npart_list;

    lennard_jones_linear();
    void test();
};

template <typename float_type>
void lennard_jones_linear<float_type>::test()
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
    BOOST_CHECK( get_en_pot(*force, en_pot.begin()) == en_pot.end() );

    std::vector<vector_type> f_list(particle->nparticle());
    BOOST_CHECK( get_net_force(*force, f_list.begin()) == f_list.end() );

    // FIXME find better estimate for floating point error, a global estimate
    // seems unsuitable due to the subtractions in the potential computation
    const float_type tolerance = 5e2 * 10 * std::numeric_limits<float_type>::epsilon();

    for (unsigned int i = 0; i < npart; ++i) {
        unsigned int type1 = species[i];
        unsigned int type2 = species[(i + 1) % npart];
        vector_type r = r_list[i] - r_list[(i + 1) % npart];
        vector_type f = f_list[i];

        // reference values from host module
        float_type fval, en_pot_;
        boost::tie(fval, en_pot_) = (*host_potential)(inner_prod(r, r), type1, type2);
        // the GPU force module stores only a fraction of these values
        en_pot_ /= 2;

        BOOST_CHECK_SMALL(norm_inf(fval * r - f), norm_inf(f) * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 4 * tolerance);
    }
}

template <typename float_type>
lennard_jones_linear<float_type>::lennard_jones_linear()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart_list = {1000, 2};
    float box_length = 100;
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }
    float cutoff = box_length / 2;

    typedef typename potential_type::matrix_type matrix_type;
    matrix_type cutoff_array(2, 2);
    matrix_type epsilon_array(2, 2);
    matrix_type sigma_array(2, 2);
    cutoff_array <<= cutoff, cutoff, cutoff, cutoff;
    epsilon_array <<= 1.f, .5f, .5f, .25f;
    sigma_array <<= 1.f, 2.f, 2.f, 4.f;

    // create modules
    particle = std::make_shared<particle_type>(
        accumulate(npart_list.begin(), npart_list.end(), 0)
      , npart_list.size()
    );
    box = std::make_shared<box_type>(edges);

    potential = std::make_shared<potential_type>(
        cutoff_array, epsilon_array, sigma_array
    );
    host_potential = std::make_shared<host_potential_type>(
        cutoff_array, epsilon_array, sigma_array
    );
    neighbour = std::make_shared<neighbour_type>(particle);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour);
}

BOOST_FIXTURE_TEST_CASE( lennard_jones_linear_gpu, halmd::device ) {
    lennard_jones_linear<float>().test();
}
#endif // HALMD_WITH_GPU
