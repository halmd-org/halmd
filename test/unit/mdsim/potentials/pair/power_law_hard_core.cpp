/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Michael Kopp
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

#define BOOST_TEST_MODULE power_law_hard_core
#include <boost/test/unit_test.hpp>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <boost/numeric/ublas/banded.hpp>
#include <cmath> // std::pow
#include <limits>
#include <numeric> // std::accumulate

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/pair/power_law_hard_core.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/shifted.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/pair/power_law_hard_core.hpp>
# include <halmd/mdsim/gpu/potentials/pair/truncations/shifted.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/pair/gpu/neighbour_chain.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/** test power law potential with divergence at a finite distance (hard core)
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

BOOST_AUTO_TEST_CASE( power_law_hard_core_host )
{
    typedef mdsim::host::potentials::pair::power_law<double> base_potential_type;
    typedef mdsim::host::potentials::pair::adapters::hard_core<base_potential_type> modified_potential_type;
    typedef mdsim::host::potentials::pair::truncations::shifted<modified_potential_type> potential_type;
    typedef potential_type::matrix_type matrix_type;
    typedef potential_type::uint_matrix_type uint_matrix_type;

    // define interaction parameters
    //
    // choose numbers that are exactly representable as float,
    // otherwise one has to account for rounding errors in the
    // computation of reference values
    unsigned int ntype = 2;  // test a binary mixture
    matrix_type cutoff_array(ntype, ntype);
    cutoff_array <<=
        5., 6.
      , 6., 7.;
    matrix_type core_array(ntype, ntype);
    core_array <<=
        0.375, 0.5
      , 0.5, 0.75;
    matrix_type epsilon_array(ntype, ntype);
    epsilon_array <<=
        1., .5
      , .5, .25;
    matrix_type sigma_array(ntype, ntype);
    sigma_array <<=
        1., 2.
      , 2., 4.;
    uint_matrix_type index_array(ntype, ntype);
    index_array <<=
        6, 12
      , 12, 7;

    // construct module
    potential_type potential(cutoff_array, core_array, epsilon_array, sigma_array, index_array);

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

    matrix_type cutoff = potential.r_cut_sigma();
    BOOST_CHECK(cutoff(0, 0) == cutoff_array(0, 0));
    BOOST_CHECK(cutoff(0, 1) == cutoff_array(0, 1));
    BOOST_CHECK(cutoff(1, 0) == cutoff_array(1, 0));
    BOOST_CHECK(cutoff(1, 1) == cutoff_array(1, 1));

    matrix_type core = potential.r_core_sigma();
    BOOST_CHECK(core(0, 0) == core_array(0, 0));
    BOOST_CHECK(core(0, 1) == core_array(0, 1));
    BOOST_CHECK(core(1, 0) == core_array(1, 0));
    BOOST_CHECK(core(1, 1) == core_array(1, 1));

    uint_matrix_type index = potential.index();
    BOOST_CHECK(index(0, 0) == index_array(0, 0));
    BOOST_CHECK(index(0, 1) == index_array(0, 1));
    BOOST_CHECK(index(1, 0) == index_array(1, 0));
    BOOST_CHECK(index(1, 1) == index_array(1, 1));

    // evaluate some points of potential and force
    typedef boost::array<double, 3> array_type;
    const double eps = numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot)
    // interaction AA: ε=1, σ=1, rc=5σ, r_core=0.375σ, n=6
    boost::array<array_type, 5> results_aa = {{
        {{0.5, 2.5165824e7, 262143.9998978285}}
      , {{0.6666666666666667, 50122.75353685236, 1624.34839207839}}
      , {{1., 161.0612736, 16.77711382854513}}
      , {{2., 0.1002646166123735, 0.05420782921016794}}
      , {{10., 7.840541955650863e-8, -0.0001009137012623106}}
    }};

    BOOST_FOREACH (array_type const& a, results_aa) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 0);  // interaction AA

        // tolerance due to floating-point rounding depends on difference (r-r_core)
        double r = a[0] / sigma_array(0, 0);        //< r in units of σ
        double tolerance = eps * index_array(0, 0) * (1 + r / (r - core_array(0, 0)));

        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction AB: ε=0.5, σ=2, rc=6σ, r_core=0.5σ, n=12
    boost::array<array_type, 5> results_ab = {{
        {{1.333333333333333, 2.9386561536e10, 1.088391168e9}}
      , {{1.5, 1.34217728e8, 8.388607999999999e6}}
      , {{2., 12288., 2047.999999999347}}
      , {{4., 0.00385367331462947, 0.003853672662073555}}
      , {{20., 2.92202811508691e-14, -6.516306057676999e-10}}
    }};

    BOOST_FOREACH (array_type const& a, results_ab) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 0, 1);  // interaction AB

        // tolerance due to floating-point rounding depends on difference (r-r_core)
        double r = a[0] / sigma_array(0, 1);        //< r in units of σ
        double tolerance = eps * index_array(0, 1) * (1 + r / (r - core_array(0, 1)));

        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };

    // interaction BB: ε=0.25, σ=4, rc=7σ, r_core=.75σ, n=7
    boost::array<array_type, 5> results_bb = {{
        {{3.333333333333333, 5.64350976e7, 8.957951999999329e6}}
      , {{3.5, 2.097152e6, 524287.9999993289}}
      , {{4., 7168., 4095.999999328911}}
      , {{8., 0.00917504, 0.05242812891136}}
      , {{20., 2.055117329014365e-7, 9.310909815212629e-6}}
    }};

    BOOST_FOREACH (array_type const& a, results_bb) {
        double rr = std::pow(a[0], 2);
        double fval, en_pot;
        std::tie(fval, en_pot) = potential(rr, 1, 1);  // interaction BB

        // tolerance due to floating-point rounding depends on difference (r-r_core)
        double r = a[0] / sigma_array(1, 1);        //< r in units of σ
        double tolerance = eps * index_array(1, 1) * (1 + r / (r - core_array(1, 1)));

        BOOST_CHECK_CLOSE_FRACTION(fval, a[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, a[2], tolerance);
    };
}

#ifdef HALMD_WITH_GPU

template <typename float_type>
struct power_law_hard_core
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::pair::power_law<float_type> base_potential_type;
    typedef mdsim::gpu::potentials::pair::adapters::hard_core<base_potential_type> modified_potential_type;
    typedef mdsim::gpu::potentials::pair::truncations::shifted<modified_potential_type> potential_type;

    typedef mdsim::host::potentials::pair::power_law<double> base_host_potential_type;
    typedef mdsim::host::potentials::pair::adapters::hard_core<base_host_potential_type> modified_host_potential_type;
    typedef mdsim::host::potentials::pair::truncations::shifted<modified_host_potential_type> host_potential_type;
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

    power_law_hard_core();
    void test();
};

template <typename float_type>
void power_law_hard_core<float_type>::test()
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

    for (unsigned int i = 0; i < npart; i += 100) {
        unsigned int type1 = species[i];
        unsigned int type2 = species[(i + 1) % npart];
        vector_type r = r_list[i] - r_list[(i + 1) % npart];
        vector_type f = f_list[i];

        // reference values from host module
        float_type fval, en_pot_;
        std::tie(fval, en_pot_) = (*host_potential)(inner_prod(r, r), type1, type2);
        // the GPU force module stores only a fraction of these values
        en_pot_ /= 2;

        // estimated upper bound on floating-point error for (r - r_core)^n
        float_type const eps = numeric_limits<float_type>::epsilon();
        unsigned int index = host_potential->index()(type1, type2);
        float_type r_core = host_potential->r_core_sigma()(type1, type2) * host_potential->sigma()(type1, type2);
        float_type r_norm = norm_2(r);
        float_type tolerance = eps * index * (1 + r_norm / (r_norm - r_core));

        // check both absolute and relative error
        BOOST_CHECK_SMALL(norm_inf(fval * r - f), max(norm_inf(fval * r), float_type(1)) * tolerance);

        // the prefactor is not justified, it is needed for absolute differences on the order 1e-14
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 3 * tolerance);
    }
}

template <typename float_type>
power_law_hard_core<float_type>::power_law_hard_core()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart_list.push_back(100);
    npart_list.push_back(2); // particle distance should be large than hard core
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
    matrix_type core_array(2, 2);
    core_array <<=
        0.375, 0.5
      , 0.5, 0.75;
    matrix_type epsilon_array(2, 2);
    epsilon_array <<=
        1., .5
      , .5, .25;
    matrix_type sigma_array(2, 2);
    sigma_array <<=
        1., 2.
      , 2., 4.;
    uint_matrix_type index_array(2, 2);
    index_array <<=
        6, 12
      , 12, 7;

    // create modules
    particle = std::make_shared<particle_type>(accumulate(npart_list.begin(), npart_list.end(), 0), npart_list.size());
    box = std::make_shared<box_type>(edges);
    potential = std::make_shared<potential_type>(
        cutoff_array, core_array, epsilon_array, sigma_array, index_array
    );
    host_potential = std::make_shared<host_potential_type>(
        cutoff_array, core_array, epsilon_array, sigma_array, index_array
    );
    neighbour = std::make_shared<neighbour_type>(particle);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour);
    particle->on_prepend_force([=](){force->check_cache();});
    particle->on_force([=](){force->apply();});
}

BOOST_FIXTURE_TEST_CASE( power_law_hard_core_gpu, device ) {
    power_law_hard_core<float>().test();
}
#endif // HALMD_WITH_GPU
