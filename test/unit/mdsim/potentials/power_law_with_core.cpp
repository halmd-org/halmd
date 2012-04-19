/*
 * Copyright © 2011-2012  Michael Kopp and Felix Höfling
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

#define BOOST_TEST_MODULE power_law_with_core
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <cmath> // std::pow
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/power_law_with_core.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_kernel.cuh>
# include <halmd/mdsim/gpu/potentials/power_law_with_core.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/gpu/neighbour_chain.hpp>
#endif


using namespace boost;
using namespace boost::assign; // list_of
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

BOOST_AUTO_TEST_CASE( power_law_with_core_host )
{
    typedef mdsim::host::potentials::power_law_with_core<double> potential_type;
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
    potential_type potential(ntype, ntype, cutoff_array, core_array, epsilon_array, sigma_array, index_array);

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

    // evaluate some points of the potential, force, and hypervirial
    typedef boost::tuple<double, double, double, double> tuple_type;
    const double eps = numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot, hvir)
    // interaction AA: ε=1, σ=1, rc=5σ, r_core=0.375σ, n=6
    vector<tuple_type> results = tuple_list_of
        (0.5, 2.5165824e7, 262143.9998978285, 1.69869312e8)
        (0.6666666666666667, 50122.75353685236, 1624.34839207839, 334151.6902456824)
        (1., 161.0612736, 16.77711382854513, 1642.82499072)
        (2., 0.1002646166123735, 0.05420782921016794, 3.054214475269223)
        (10., 7.840541955650863e-8, -0.0001009137012623106, 0.0000491815813581736);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 0);  // interaction AA

        // tolerance due to floating-point rounding depends on difference (r-r_core)
        double r = get<0>(a) / sigma_array(0, 0);        //< r in units of σ
        double tolerance = eps * index_array(0, 0) * (1 + r / (r - core_array(0, 0)));

        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction AB: ε=0.5, σ=2, rc=6σ, r_core=0.5σ, n=12
    results = tuple_list_of
        (1.333333333333333, 2.9386561536e10, 1.088391168e9, 2.664381579264e12)
        (1.5, 1.34217728e8, 8.388607999999999e6, 1.1475615744e10)
        (2., 12288., 2047.999999999347, 1.2288e6)
        (4., 0.00385367331462947, 0.003853672662073555, 1.007093292889835)
        (20., 2.92202811508691e-14, -6.516306057676999e-10, 1.482544791023043e-10);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 1);  // interaction AB

        // tolerance due to floating-point rounding depends on difference (r-r_core)
        double r = get<0>(a) / sigma_array(0, 1);        //< r in units of σ
        double tolerance = eps * index_array(0, 1) * (1 + r / (r - core_array(0, 1)));

        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction BB: ε=0.25, σ=4, rc=7σ, r_core=.75σ, n=7
    results = tuple_list_of
        (3.333333333333333, 5.64350976e7, 8.957951999999329e6, 4.953747456e10)
        (3.5, 2.097152e6, 524287.9999993289, 1.41295616e9)
        (4., 7168., 4095.999999328911, 3.555328e6)
        (8., 0.00917504, 0.05242812891136, 6.928990208)
        (20., 2.055117329014365e-7, 9.310909815212629e-6, 0.0006914865365860098);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 1, 1);  // interaction BB

        // tolerance due to floating-point rounding depends on difference (r-r_core)
        double r = get<0>(a) / sigma_array(1, 1);        //< r in units of σ
        double tolerance = eps * index_array(1, 1) * (1 + r / (r - core_array(1, 1)));

        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };
}

#ifdef WITH_CUDA

template <typename float_type>
struct power_law_with_core
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::power_law_with_core<float_type> potential_type;
    typedef mdsim::host::potentials::power_law_with_core<double> host_potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef test::unit::mdsim::potentials::gpu::neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    shared_ptr<box_type> box;
    shared_ptr<potential_type> potential;
    shared_ptr<force_type> force;
    shared_ptr<neighbour_type> neighbour;
    shared_ptr<particle_type> particle;
    shared_ptr<host_potential_type> host_potential;
    vector<unsigned int> npart_list;

    power_law_with_core();
    void test();
};

template <typename float_type>
void power_law_with_core<float_type>::test()
{
    using namespace halmd::mdsim::gpu; // particle_kernel::{tagged, untagged}

    // place particles along the x-axis within one half of the box,
    // put every second particle at the origin
    unsigned int npart = particle->nparticle();
    vector_type dx(0);
    dx[0] = box->edges()[0][0] / npart / 2;

    cuda::host::vector<float4> r_list(particle->position().size());
    for (unsigned int k = 0; k < r_list.size(); ++k) {
        vector_type r = (k % 2) ? k * dx : vector_type(0);
        unsigned int type = (k < npart_list[0]) ? 0U : 1U;  // set particle type for a binary mixture
        r_list[k] = particle_kernel::tagged<vector_type>(r, type);
    }
    cuda::copy(r_list, particle->position());

    particle->aux_enable();              // enable computation of auxiliary quantities
    particle->prepare();
    force->compute();

    // read forces and other stuff from device
    cuda::host::vector<typename particle_type::gpu_vector_type> f_list(particle->force().size());
    cuda::copy(particle->force(), f_list);

    cuda::host::vector<float> en_pot(particle->en_pot().size());
    cuda::copy(particle->en_pot(), en_pot);

    cuda::host::vector<float> hypervirial(particle->hypervirial().size());
    cuda::copy(particle->hypervirial(), hypervirial);

    for (unsigned int i = 0; i < npart; i += 100) {
        vector_type r1, r2;
        unsigned int type1, type2;
        tie(r1, type1) = particle_kernel::untagged<vector_type>(r_list[i]);
        tie(r2, type2) = particle_kernel::untagged<vector_type>(r_list[(i + 1) % npart]);
        vector_type r = r1 - r2;
        vector_type f = f_list[i];

        // reference values from host module
        float_type fval, en_pot_, hvir;
        tie(fval, en_pot_, hvir) = (*host_potential)(inner_prod(r, r), type1, type2);
        // the GPU force module stores only a fraction of these values
        en_pot_ /= 2;
        hvir /= 2 * dimension * dimension;

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
        BOOST_CHECK_CLOSE_FRACTION(hvir, hypervirial[i], tolerance);
    }
}

template <typename float_type>
power_law_with_core<float_type>::power_law_with_core()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart_list = list_of(100)(2); // particle distance should be large than hard core
    vector<double> mass = list_of(1)(1);
    float box_length = 100;
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
    particle = make_shared<particle_type>(npart_list, mass);
    box = make_shared<box_type>(typename box_type::vector_type(box_length));
    potential = make_shared<potential_type>(
        particle->nspecies(), particle->nspecies(), cutoff_array, core_array, epsilon_array, sigma_array, index_array
    );
    host_potential = make_shared<host_potential_type>(
        particle->nspecies(), particle->nspecies(), cutoff_array, core_array, epsilon_array, sigma_array, index_array
    );
    neighbour = make_shared<neighbour_type>(particle);
    force = make_shared<force_type>(potential, particle, particle, box, neighbour);
}

BOOST_FIXTURE_TEST_CASE( power_law_with_core_gpu, device ) {
    power_law_with_core<float>().test();
}
#endif // WITH_CUDA
