/*
 * Copyright © 2011-2023 Felix Höfling
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE morse
#include <boost/test/unit_test.hpp>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <boost/numeric/ublas/banded.hpp>
#include <cmath> // std::pow
#include <limits>
#include <numeric> // std::accumulate

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/morse.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/shifted.hpp>
# include <halmd/mdsim/gpu/potentials/pair/morse.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/pair/gpu/neighbour_chain.hpp>
# include <test/tools/cuda.hpp>
#endif
#include <test/tools/ctest.hpp>
#include <test/tools/dsfloat.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/** test the (truncated and shifted) Morse potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

BOOST_AUTO_TEST_CASE( morse_host )
{
#ifndef USE_HOST_SINGLE_PRECISION
    typedef double float_type;
#else
    typedef float float_type;
#endif
    typedef mdsim::host::potentials::pair::morse<float_type> base_potential_type;
    typedef mdsim::host::potentials::pair::truncations::shifted<base_potential_type> potential_type;
    typedef potential_type::matrix_type matrix_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    matrix_type cutoff_array(ntype, ntype);
    cutoff_array <<=
        5., 5.
      , 5., 8.;
    matrix_type epsilon_array(ntype, ntype);
    epsilon_array <<=
        1., .5
      , .5, .5;
    matrix_type sigma_array(ntype, ntype);
    sigma_array <<=
        1., 2.
      , 2., .75;
    matrix_type r_min_array(ntype, ntype);
    r_min_array <<=
        .25, 1.5
      , 1.5, 2.;
    matrix_type distortion_array(ntype, ntype);
    distortion_array <<=
        1.41, 1.
      , 1., 1.5;

    // construct module
    potential_type potential(cutoff_array, epsilon_array, sigma_array, r_min_array, distortion_array);

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

    matrix_type r_min_sigma = potential.r_min_sigma();
    BOOST_CHECK(r_min_sigma(0, 0) == r_min_array(0, 0) / sigma_array(0, 0));
    BOOST_CHECK(r_min_sigma(0, 1) == r_min_array(0, 1) / sigma_array(0, 1));
    BOOST_CHECK(r_min_sigma(1, 0) == r_min_array(1, 0) / sigma_array(1, 0));
    BOOST_CHECK(r_min_sigma(1, 1) == r_min_array(1, 1) / sigma_array(1, 1));

    matrix_type distortion = potential.distortion();
    BOOST_CHECK(distortion(0, 0) == distortion_array(0, 0));
    BOOST_CHECK(distortion(0, 1) == distortion_array(0, 1));
    BOOST_CHECK(distortion(1, 0) == distortion_array(1, 0));
    BOOST_CHECK(distortion(1, 1) == distortion_array(1, 1));

    // evaluate some points of potential and force
    typedef boost::array<float_type, 3> array_type;
    float_type const eps = numeric_limits<float_type>::epsilon();

    // expected results (r, fval, en_pot) for ε=1, σ=1, rmin=.25, B=1.41, rc=5σ
    boost::array<array_type, 6> results_aa = {{
        {{0.25, 0.0, -0.9540005656159722}}
      , {{0.5, -0.6507845119203753, -0.9069122207759784}}
      , {{0.75, -0.577737379921433, -0.8091011280363248}}
      , {{1.0, -0.4423441678610893, -0.6983391318532663}}
      , {{2.0, -0.1335362022766309, -0.3377631169819678}}
      , {{5.0, -0.006524526474499887, 0.0}}
    }};

    BOOST_FOREACH (array_type const& a, results_aa) {
        float_type rr = std::pow(a[0], 2);
        float_type fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 0);  // interaction AA

        float_type tolerance = 2 * eps;
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], 2 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction AB: ε=.5, σ=2, rmin=1.5, B=1, rc=5σ
    boost::array<array_type, 7> results_ab = {{
        {{0.25, 3.244194000059237, -0.1089119789768079}}
      , {{0.5, 1.069560557758917, -0.2754178567461116}}
      , {{0.75, 0.4413390679963154, -0.3823289065873701}}
      , {{1.0, 0.1823479270061933, -0.4455022816131835}}
      , {{2.0, -0.04306753083969286, -0.4613729534905943}}
      , {{5.0, -0.01435765600281266, -0.144512752014792}}
      , {{10.0, -0.0007030382769994306, 0.0}}
    }};

    BOOST_FOREACH (array_type const& a, results_ab) {
        float_type rr = std::pow(a[0], 2);
        float_type fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 1);  // interaction AB

        float_type tolerance = 2.5 * eps;
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], 2 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction BB: ε=.5, σ=0.75, rmin=2, B=1.5, rc=8σ
    boost::array<array_type, 7> results_bb = {{
        {{0.5, 456.725857826345, 55.21226180109373}}
      , {{0.75, 110.7622303429022, 19.2674158206776}}
      , {{1.0, 29.80898547363259, 6.254397161406729}}
      , {{2.0, 0.0, -0.481636479857814}}
      , {{3.0, -0.07481840983208633, -0.2433064324854725}}
      , {{5.0, -0.007940263658336938, -0.02630353504213594}}
      , {{6.0, -0.002720513166607205, 0.0}}
    }};

    BOOST_FOREACH (array_type const& a, results_bb) {
        float_type rr = std::pow(a[0], 2);
        float_type fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 1, 1);  // interaction BB

        float_type tolerance = 4 * eps;
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], 2 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };
}

#ifdef HALMD_WITH_GPU

template <typename float_type>
struct morse
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::pair::morse<float> base_potential_type;
    typedef mdsim::gpu::potentials::pair::truncations::shifted<base_potential_type> potential_type;
#ifndef USE_HOST_SINGLE_PRECISION
    typedef mdsim::host::potentials::pair::morse<double> base_host_potential_type;
#else
    typedef mdsim::host::potentials::pair::morse<float> base_host_potential_type;
#endif
    typedef mdsim::host::potentials::pair::truncations::shifted<base_host_potential_type> host_potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    std::shared_ptr<box_type> box;
    std::shared_ptr<potential_type> potential;
    std::shared_ptr<force_type> force;
    std::shared_ptr<neighbour_type> neighbour;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<host_potential_type> host_potential;
    vector<unsigned int> npart_list;

    morse();
    void test();
};

template <typename float_type>
void morse<float_type>::test()
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

    for (unsigned int i = 0; i < npart; ++i) {
        unsigned int type1 = species[i];
        unsigned int type2 = species[(i + 1) % npart];
        vector_type r = r_list[i] - r_list[(i + 1) % npart];
        vector_type f = f_list[i];

        // reference values from host module
        float_type fval(0), en_pot_(0);
        if (host_potential->within_range(inner_prod(r, r), type1, type2)) {
            std::tie(fval, en_pot_) = (*host_potential)(inner_prod(r, r), type1, type2);
            // the GPU force module stores only a fraction of these values
            en_pot_ /= 2;
        }

        // rough upper bound on floating-point error
        float_type const eps = numeric_limits<float>::epsilon();
        float_type tolerance = 2.5 * eps;

        // check both absolute and relative error
        BOOST_CHECK_SMALL(float_type(norm_inf(fval * r - f)), max(float_type(norm_inf(fval * r)), float_type(1)) * tolerance);

        BOOST_CHECK_SMALL(abs(double(en_pot_ - en_pot[i])), max(abs(double(en_pot_)), 1.) * double(tolerance));
    }
}

template <typename float_type>
morse<float_type>::morse()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart_list.push_back(1000);
    npart_list.push_back(2);
    float box_length = 50;
    unsigned int const dimension = box_type::vector_type::static_size;
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }
    float cutoff = box_length / 2;

    typedef typename potential_type::matrix_type matrix_type;
    matrix_type cutoff_array(2, 2);
    cutoff_array <<=
        cutoff, cutoff
      , cutoff, cutoff;
    matrix_type epsilon_array(2, 2);
    epsilon_array <<=
        1., .5
      , .5, .25;
    matrix_type sigma_array(2, 2);
    sigma_array <<=
        1., 2.
      , 2., 4.;
    matrix_type r_min_array(2, 2);
    r_min_array <<=
         .5, 1.5
      , 1.5, 1.;
    matrix_type distortion_array(2, 2);
    distortion_array <<=
        1.,  1.5
      , 1.5, 2.4;

    // create modules
    particle = std::make_shared<particle_type>(accumulate(npart_list.begin(), npart_list.end(), 0), npart_list.size());
    box = std::make_shared<box_type>(edges);
    potential = std::make_shared<potential_type>(
        cutoff_array, epsilon_array, sigma_array, r_min_array, distortion_array
    );
    host_potential = std::make_shared<host_potential_type>(
        cutoff_array, epsilon_array, sigma_array, r_min_array, distortion_array
    );
    neighbour = std::make_shared<neighbour_type>(particle);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour);
    particle->on_prepend_force([=](){force->check_cache();});
    particle->on_force([=](){force->apply();});
}

# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( morse_gpu_dsfloat, set_cuda_device ) {
    morse<dsfloat>().test();
}
# endif
# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( morse_gpu_float, set_cuda_device ) {
    morse<float>().test();
}
# endif
#endif // HALMD_WITH_GPU
