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

#define BOOST_TEST_MODULE mie
#include <boost/test/unit_test.hpp>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <boost/numeric/ublas/banded.hpp>
#include <cmath> // std::pow
#include <limits>
#include <numeric> // std::accumulate

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/mie.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/shifted.hpp>
# include <halmd/mdsim/gpu/potentials/pair/mie.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/pair/gpu/neighbour_chain.hpp>
# include <test/tools/cuda.hpp>
#endif
#include <test/tools/ctest.hpp>
#include <test/tools/dsfloat.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/** test the (truncated and shifted) Mie potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

BOOST_AUTO_TEST_CASE( mie_host )
{
#ifndef USE_HOST_SINGLE_PRECISION
    typedef double float_type;
#else
    typedef float float_type;
#endif
    typedef mdsim::host::potentials::pair::mie<float_type> base_potential_type;
    typedef mdsim::host::potentials::pair::truncations::shifted<base_potential_type> potential_type;
    typedef potential_type::matrix_type matrix_type;
    typedef potential_type::uint_matrix_type uint_matrix_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    matrix_type cutoff_array(ntype, ntype);
    cutoff_array <<=
        5., 5.
      , 5., 5.;
    matrix_type epsilon_array(ntype, ntype);
    epsilon_array <<=
        1., .5
      , .5, .25;
    matrix_type sigma_array(ntype, ntype);
    sigma_array <<=
        1., 2.
      , 2., 4.;
    uint_matrix_type index_m_array(ntype, ntype);
    index_m_array <<=
        12, 4
      ,  4, 12;
    uint_matrix_type index_n_array(ntype, ntype);
    index_n_array <<=
        2, 2
      , 2, 6;

    // construct module
    potential_type potential(cutoff_array, epsilon_array, sigma_array, index_m_array, index_n_array);

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

    uint_matrix_type index_m = potential.index_m();
    BOOST_CHECK(index_m(0, 0) == index_m_array(0, 0));
    BOOST_CHECK(index_m(0, 1) == index_m_array(0, 1));
    BOOST_CHECK(index_m(1, 0) == index_m_array(1, 0));
    BOOST_CHECK(index_m(1, 1) == index_m_array(1, 1));

    uint_matrix_type index_n = potential.index_n();
    BOOST_CHECK(index_n(0, 0) == index_n_array(0, 0));
    BOOST_CHECK(index_n(0, 1) == index_n_array(0, 1));
    BOOST_CHECK(index_n(1, 0) == index_n_array(1, 0));
    BOOST_CHECK(index_n(1, 1) == index_n_array(1, 1));

    // evaluate some points of potential and force
    typedef boost::array<float_type, 3> array_type;
    const float_type tolerance = 6 * numeric_limits<float_type>::epsilon();

    // expected results (r, fval, en_pot) for ε=1, σ=1, m=12, n=2, rc=5σ
    boost::array<array_type, 5> results_aa = {{
        {{0.5, 337553.0137048161, 7026.699262368106}}
      , {{0.8, 460.1405827001969, 22.37362839697643}}
      , {{1.0, 17.17162897326307, 0.06868650885955306}}
      , {{2.0, -0.2133876744968482, -0.3601849862490435}}
      , {{10.0, -0.0003434325792592018, 0.05151487988800715}}
    }};

    BOOST_FOREACH (array_type const& a, results_aa) {
        float_type rr = std::pow(a[0], 2);
        float_type fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 0);  // interaction AA
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction AB: ε=.5, σ=2, m=4, n=2, rc=5σ
    boost::array<array_type, 5> results_ab = {{
        {{0.5, 7936.0, 480.0768}}
      , {{0.8, 449.21875, 65.7018}}
      , {{1.0, 112.0, 24.0768}}
      , {{2.0, 1.0, 0.0768}}
      , {{10.0, -0.001472, 0.0}}
    }};

    BOOST_FOREACH (array_type const& a, results_ab) {
        float_type rr = std::pow(a[0], 2);
        float_type fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 1);  // interaction AB
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction BB: ε=.25, σ=4, m=12, n=6, rc=5σ
    boost::array<array_type, 5> results_bb = {{
        {{0.5, 3298528591872.0, 68719214592.00006}}
      , {{0.8, 4577490234.375, 244125000.000064}}
      , {{1.0, 201302016.0, 16773120.000064}}
      , {{2.0, 12192.0, 4032.000063995904}}
      , {{10.0, -0.00024374673408, -0.00401522688}}
    }};

    BOOST_FOREACH (array_type const& a, results_bb) {
        float_type rr = std::pow(a[0], 2);
        float_type fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 1, 1);  // interaction BB
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };
}

#ifdef HALMD_WITH_GPU

template <typename float_type>
struct mie
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::pair::mie<float> base_potential_type;
    typedef mdsim::gpu::potentials::pair::truncations::shifted<base_potential_type> potential_type;
#ifndef USE_HOST_SINGLE_PRECISION
    typedef mdsim::host::potentials::pair::mie<double> base_host_potential_type;
#else
    typedef mdsim::host::potentials::pair::mie<float> base_host_potential_type;
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

    mie();
    void test();
};

template <typename float_type>
void mie<float_type>::test()
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

    const float_type tolerance = 20 * numeric_limits<float>::epsilon(); // FIXME the prefactor is an unjustified guess

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

        // FIXME the tolerance needs to cover both very large and vanishing forces
        BOOST_CHECK_SMALL(float_type (norm_inf(fval * r - f)), max(norm_inf(fval * r), float_type(1)) * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 4 * tolerance);
    }
}

template <typename float_type>
mie<float_type>::mie()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart_list.push_back(1000);
    npart_list.push_back(2);
    float box_length = 100;
    unsigned int const dimension = box_type::vector_type::static_size;
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }
    float cutoff = box_length / 2;

    typedef typename potential_type::matrix_type matrix_type;
    typedef typename potential_type::uint_matrix_type uint_matrix_type;
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
    uint_matrix_type index_m_array(2, 2);
    index_m_array <<=
        12, 12
      , 12, 12;
    uint_matrix_type index_n_array(2, 2);
    index_n_array <<=
        4, 2
      , 2, 6;

    // create modules
    particle = std::make_shared<particle_type>(accumulate(npart_list.begin(), npart_list.end(), 0), npart_list.size());
    box = std::make_shared<box_type>(edges);
    potential = std::make_shared<potential_type>(
        cutoff_array, epsilon_array, sigma_array, index_m_array, index_n_array
    );
    host_potential = std::make_shared<host_potential_type>(
        cutoff_array, epsilon_array, sigma_array, index_m_array, index_n_array
    );
    neighbour = std::make_shared<neighbour_type>(particle);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour);
    particle->on_prepend_force([=](){force->check_cache();});
    particle->on_force([=](){force->apply();});
}

# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( mie_gpu_dsfloat, set_cuda_device ) {
    mie<dsfloat>().test();
}
# endif
# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( mie_gpu_float, set_cuda_device ) {
    mie<float>().test();
}
# endif
#endif // HALMD_WITH_GPU
