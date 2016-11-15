/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2016      Daniel Kirchner
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

#define BOOST_TEST_MODULE lennard_jones
#include <boost/test/unit_test.hpp>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <boost/numeric/ublas/banded.hpp>
#include <cmath> // std::pow
#include <limits>
#include <numeric> // std::accumulate

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/force_shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/sharp.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/smooth_r4.hpp>
#include <halmd/mdsim/host/potentials/pair/lennard_jones.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/force_shifted.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/sharp.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/shifted.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/smooth_r4.hpp>
# include <halmd/mdsim/gpu/potentials/pair/lennard_jones.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/pair/gpu/neighbour_chain.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/** test Lennard-Jones potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

/**
 * expected results for the host module
 *
 * This template is specialized for every available potential truncation and
 * provides the expected results for the potential test.
 */
template<typename T>
struct results;

template<>
struct results<mdsim::host::potentials::pair::truncations::sharp<mdsim::host::potentials::pair::lennard_jones<double>>>
{
    typedef boost::array<double, 3> array_type;

    static constexpr boost::array<array_type, 5> aa() {
        return {{
            {{0.2, 2.92959375e11, 9.7649999999999905e8}}
          , {{0.5, 780288., 16128.}}
          , {{1., 24., 0.}}
          , {{2., -0.0908203125, -0.0615234375}}
          , {{10., -2.3999952e-7, -3.9999960000000006e-06}}
        }};
    }
    static constexpr boost::array<array_type, 5> ab() {
        return {{
            {{0.2, 5.999997e14, 1.9999980000000002e12}}
          , {{0.5, 1.610416128e9, 3.354624e7}}
          , {{1., 97536., 8064.}}
          , {{2., 3., 0.}}
          , {{10., -7.67901696e-6,-0.00012799180800000004}}
        }};
    }
    static constexpr boost::array<array_type, 5> bb() {
        return {{
            {{0.2, 1.2287999904e18, 4.095999936e15}}
          , {{0.5, 3.298528591872e12, 6.871921459200006e10}}
          , {{1., 2.01302016e8, 1.677312e7}}
          , {{2., 12192., 4032.}}
          , {{10., -0.00024374673408, -0.0040792227840000007}}
        }};
    }
};

template<>
struct results<mdsim::host::potentials::pair::truncations::shifted<mdsim::host::potentials::pair::lennard_jones<double>>>
{
    typedef boost::array<double, 3> array_type;

    static constexpr boost::array<array_type, 5> aa() {
        return {{
            {{0.2, 2.92959375e11, 9.76500000000256e8}}
          , {{0.5, 780288., 16128.00025598362}}
          , {{1., 24., 0.000255983616}}
          , {{2., -0.0908203125, -0.061267453884}}
          , {{10., -2.3999952e-7, 0.00025198362}}
        }};
    }
    static constexpr boost::array<array_type, 5> ab() {
        return {{
            {{0.2, 5.999997e14, 1.9999980000000002e12}}
          , {{0.5, 1.610416128e9, 3.3546240000127994e7}}
          , {{1., 97536., 8064.000127991808}}
          , {{2., 3., 0.000127991808}}
          , {{10., -7.67901696e-6,0.}}
        }};
    }
    static constexpr boost::array<array_type, 5> bb() {
        return {{
            {{0.2, 1.2287999904e18, 4.095999936e15}}
          , {{0.5, 3.298528591872e12, 6.871921459200006e10}}
          , {{1., 2.01302016e8, 1.6773120000064e7}}
          , {{2., 12192., 4032.000063995904}}
          , {{10., -0.00024374673408, -0.00401522688}}
        }};
    }
};

template<>
struct results<mdsim::host::potentials::pair::truncations::force_shifted<mdsim::host::potentials::pair::lennard_jones<double>>>
{
    typedef boost::array<double, 3> array_type;

    static constexpr boost::array<array_type, 5> aa() {
        return {{
            {{0.2, 2.929593750000015e11, 9.765000000017304e8}}
          , {{0.5, 780288.0006143214, 16128.00163820667}}
          , {{1., 24.0003071606784, 0.0014846263296}}
          , {{2., -0.0906667321608, -0.0603459718488}}
          , {{10., 0.00003047606832, -0.001283819772}}
        }};
    }
    static constexpr boost::array<array_type, 5> ab() {
        return {{
            {{0.2, 5.999997e14, 1.999998000000001e12}}
          , {{0.5, 1.610416128000154e9, 3.35462400008575e7}}
          , {{1., 97536.00007679017, 8064.000819103334}}
          , {{2., 3.0000383950848, 0.0007423131648}}
          , {{10., 0, 0}}
        }};
    }
    static constexpr boost::array<array_type, 5> bb() {
        return {{
            {{0.2, 1.2287999904e18, 4.095999936e15}}
          , {{0.5, 3.298528591872e12, 6.871921459200044e10}}
          , {{1., 2.013020160000192e8, 1.677312000042875e7}}
          , {{2., 12192.00000959877, 4032.000409551667}}
          , {{10., -0.00024182697984, -0.003823251456}}
        }};
    }
};

template<>
struct results<mdsim::host::potentials::pair::truncations::smooth_r4<mdsim::host::potentials::pair::lennard_jones<double>>>
{
    typedef boost::array<double, 3> array_type;

    static constexpr boost::array<array_type, 5> aa() {
        return {{
            {{0.2, 2.9295937499965955e11, 9.7649999999910533e8}}
          , {{0.5, 780287.99999885436, 16128.000255959034}}
          , {{1., 23.999999999941405, 0.00025598361599937509}}
          , {{2., -0.090820312499614378, -0.061267453883527251}}
          , {{10., -2.3999952001991871e-07, 0.00025198361999974805}}
        }};
    }
    static constexpr boost::array<array_type, 5> ab() {
        return {{
            {{0.2, 5.99999699999959e14, 1.9999979999998628e12}}
          , {{0.5, 1.6104161279998786e9, 3.3546240000125419e7}}
          , {{1., 97535.999999991051, 8064.0001279910393}}
          , {{2., 2.9999999999995426, 0.00012799180799998052}}
          , {{10., 0., 0.}}
        }};
    }
    static constexpr boost::array<array_type, 5> bb() {
        return {{
            {{0.2, 1.2287999903999936e18, 4.095999935999979e15}}
          , {{0.5, 3.2985285918719858e12, 6.8719214591999763e10}}
          , {{1., 2.0130201599999908e8, 1.6773120000063917e7}}
          , {{2., 12191.999999999929, 4032.0000639958798}}
          , {{10., -0.00024374673407999485, -0.0040152268799997504}}
        }};
    }
};

/**
 * To be able to create different kinds of potential trunctations, a
 * template struct make_potential is created and specialized for any
 * potential truncation type that needs additional arguments.
 */
template<typename potential_type>
struct make_potential
{
    typedef typename potential_type::matrix_type matrix_type;
    static potential_type make(matrix_type const& cutoff, matrix_type const& epsilon, matrix_type const& sigma)
    {
        return potential_type(cutoff, epsilon, sigma);
    }
};

template<typename base_potential_type>
struct make_potential<mdsim::host::potentials::pair::truncations::smooth_r4<base_potential_type>>
{
    typedef mdsim::host::potentials::pair::truncations::smooth_r4<base_potential_type> potential_type;
    typedef typename potential_type::matrix_type matrix_type;
    static potential_type make(matrix_type const& cutoff, matrix_type const& epsilon, matrix_type const& sigma)
    {
        return potential_type(cutoff, 0.005, epsilon, sigma);
    }
};

template<typename base_potential_type>
struct make_potential<mdsim::gpu::potentials::pair::truncations::smooth_r4<base_potential_type>>
{
    typedef mdsim::gpu::potentials::pair::truncations::smooth_r4<base_potential_type> potential_type;
    typedef typename potential_type::matrix_type matrix_type;
    static potential_type make(matrix_type const& cutoff, matrix_type const& epsilon, matrix_type const& sigma)
    {
        return potential_type(cutoff, 0.005, epsilon, sigma);
    }
};

/**
 * The tolerance for the GPU tests can be overridden for a specific
 * kind of potential trunctation.
 *
 * This is necessary as the subtraction in the force_shifted potential
 * truncation is ill-conditioned.
 */
template<typename float_type, typename potential_type>
struct gpu_tolerance
{
    static constexpr float_type value = 10 * numeric_limits<float_type>::epsilon();
};

template<typename float_type, typename base_potential_type>
struct gpu_tolerance<float_type, mdsim::gpu::potentials::pair::truncations::force_shifted<base_potential_type>>
{
    // FIXME find better estimate for floating point error, a global estimate
    // seems unsuitable due to the subtractions in the potential computation
    static constexpr float_type value = 5e2 * 10 * numeric_limits<float_type>::epsilon();
};

BOOST_AUTO_TEST_CASE( lennard_jones_host )
{
    typedef mdsim::host::potentials::pair::lennard_jones<double> base_potential_type;
    typedef mdsim::host::potentials::pair::truncations::TRUNCATION_TYPE<base_potential_type> potential_type;
    typedef potential_type::matrix_type matrix_type;

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

    // construct module
    potential_type potential = make_potential<potential_type>::make(cutoff_array, epsilon_array, sigma_array);

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
    typedef boost::array<double, 3> array_type;
    const double tolerance = 5 * numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot) for ε=1, σ=1, rc=5σ
    boost::array<array_type, 5> const& results_aa = results<potential_type>::aa();

    BOOST_FOREACH (array_type const& a, results_aa) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 0);  // interaction AA
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction AB: ε=.5, σ=2, rc=5σ
    boost::array<array_type, 5> const& results_ab = results<potential_type>::ab();

    BOOST_FOREACH (array_type const& a, results_ab) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 1);  // interaction AB
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction BB: ε=.25, σ=4, rc=5σ
    boost::array<array_type, 5> const& results_bb = results<potential_type>::bb();

    BOOST_FOREACH (array_type const& a, results_bb) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 1, 1);  // interaction BB
        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };
}

#ifdef HALMD_WITH_GPU

template <typename float_type>
struct lennard_jones
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::pair::lennard_jones<float_type> base_potential_type;
    typedef mdsim::gpu::potentials::pair::truncations::TRUNCATION_TYPE<base_potential_type> potential_type;
    typedef mdsim::host::potentials::pair::lennard_jones<double> base_host_potential_type;
    typedef mdsim::host::potentials::pair::truncations::TRUNCATION_TYPE<base_host_potential_type> host_potential_type;
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

    lennard_jones();
    void test();
};

template <typename float_type>
void lennard_jones<float_type>::test()
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

    const float_type tolerance = gpu_tolerance<float_type, potential_type>::value;

    for (unsigned int i = 0; i < npart; ++i) {
        unsigned int type1 = species[i];
        unsigned int type2 = species[(i + 1) % npart];
        vector_type r = r_list[i] - r_list[(i + 1) % npart];
        vector_type f = f_list[i];

        // reference values from host module
        float_type fval, en_pot_;
        std::tie(fval, en_pot_) = (*host_potential)(inner_prod(r, r), type1, type2);
        // the GPU force module stores only a fraction of these values
        en_pot_ /= 2;

        BOOST_CHECK_SMALL(norm_inf(fval * r - f), norm_inf(f) * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 4 * tolerance);
    }
}

template <typename float_type>
lennard_jones<float_type>::lennard_jones()
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
    matrix_type cutoff_array(2, 2);
    cutoff_array <<=
        cutoff, cutoff
      , cutoff, cutoff;
    matrix_type epsilon_array(2, 2);
    epsilon_array <<=
        1., .5
      , 5., .25f;
    matrix_type sigma_array(2, 2);
    sigma_array <<=
        1., 2.
      , 2., 4.;

    // create modules
    particle = std::make_shared<particle_type>(accumulate(npart_list.begin(), npart_list.end(), 0), npart_list.size());
    box = std::make_shared<box_type>(edges);
    potential = std::make_shared<potential_type>(make_potential<potential_type>::make(cutoff_array, epsilon_array, sigma_array));
    host_potential = std::make_shared<host_potential_type>(make_potential<host_potential_type>::make(cutoff_array, epsilon_array, sigma_array));
    neighbour = std::make_shared<neighbour_type>(particle);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour);
    particle->on_prepend_force([=](){force->check_cache();});
    particle->on_force([=](){force->apply();});
}

BOOST_FIXTURE_TEST_CASE( lennard_jones_gpu, device ) {
    lennard_jones<float>().test();
}
#endif // HALMD_WITH_GPU
