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

#define BOOST_TEST_MODULE modified_lennard_jones
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <cmath> // std::pow
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/modified_lennard_jones.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_kernel.cuh>
# include <halmd/mdsim/gpu/potentials/modified_lennard_jones.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/gpu/neighbour_chain.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace std;

/** test modified Lennard-Jones potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

BOOST_AUTO_TEST_CASE( modified_lennard_jones_host )
{
    typedef mdsim::host::potentials::modified_lennard_jones<double> potential_type;
    typedef potential_type::matrix_type matrix_type;
    typedef potential_type::uint_matrix_type uint_matrix_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    boost::array<float, 3> cutoff_array = list_of(5.f)(5.f)(5.f);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(4.f);
    boost::array<unsigned, 3> index_m_array = list_of(12)(12)(12);
    boost::array<unsigned, 3> index_n_array = list_of(4)(2)(6);

    // construct module
    potential_type potential(ntype, cutoff_array, epsilon_array, sigma_array, index_m_array, index_n_array);

    // test paramters
    matrix_type epsilon = potential.epsilon();
    BOOST_CHECK(epsilon(0, 0) == epsilon_array[0]);
    BOOST_CHECK(epsilon(0, 1) == epsilon_array[1]);
    BOOST_CHECK(epsilon(1, 0) == epsilon_array[1]);
    BOOST_CHECK(epsilon(1, 1) == epsilon_array[2]);

    matrix_type sigma = potential.sigma();
    BOOST_CHECK(sigma(0, 0) == sigma_array[0]);
    BOOST_CHECK(sigma(0, 1) == sigma_array[1]);
    BOOST_CHECK(sigma(1, 0) == sigma_array[1]);
    BOOST_CHECK(sigma(1, 1) == sigma_array[2]);

    uint_matrix_type index_m = potential.index_m();
    BOOST_CHECK(index_m(0, 0) == index_m_array[0]);
    BOOST_CHECK(index_m(0, 1) == index_m_array[1]);
    BOOST_CHECK(index_m(1, 0) == index_m_array[1]);
    BOOST_CHECK(index_m(1, 1) == index_m_array[2]);

    uint_matrix_type index_n = potential.index_n();
    BOOST_CHECK(index_n(0, 0) == index_n_array[0]);
    BOOST_CHECK(index_n(0, 1) == index_n_array[1]);
    BOOST_CHECK(index_n(1, 0) == index_n_array[1]);
    BOOST_CHECK(index_n(1, 1) == index_n_array[2]);

    // evaluate some points of the potential, force, and hypervirial
    typedef boost::tuple<double, double, double, double> tuple_type;
    const double tolerance = 5 * numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot, hvir) for ε=1, σ=1, m=12, n=4, rc=5σ
    vector<tuple_type> results = tuple_list_of
        (0.2, 2.929685e11, 9.765600000064e8, 1.4062496e11)
        (0.5, 785408., 16320.00639998362, 2.358272e6)
        (1., 32., 0.006399983616, 512.)
        (2., -0.2470703125, -0.242623453884, -3.859375)
        (10., -0.00001599999952, 0.00599998362, -0.006399999424);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 0);  // interaction AA
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction AB: ε=.5, σ=2, m=12, n=2, rc=5σ
    results = tuple_list_of
        (0.2, 5.9999999999e14, 1.99999999980008e12, 2.879999999992e14)
        (0.5, 1.61061248e9, 3.355440007999999e7, 4.83183808e9)
        (1., 98288., 8184.079999991808, 1.179616e6)
        (2., 5., 0.079999991808, 280.)
        (10., -0.00159999901696, 0., -0.319998820352);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 1);  // interaction AB
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction BB: ε=.25, σ=4, m=12, n=6, rc=5σ
    results = tuple_list_of
        (0.2, 1.2287999904e18, 4.095999936e15, 5.89823997696e17)
        (0.5, 3.298528591872e12, 6.871921459200006e10, 9.8955952128e12)
        (1., 2.01302016e8, 1.6773120000063997e7, 2.415771648e9)
        (2., 12192., 4032.000063995904, 587520.)
        (10., -0.00024374673408, -0.00401522688, -0.14504008089599998);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 1, 1);  // interaction BB
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };
}

#ifdef WITH_CUDA

template <typename float_type>
struct modified_lennard_jones
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::modified_lennard_jones<float_type> potential_type;
    typedef mdsim::host::potentials::modified_lennard_jones<double> host_potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef test::unit::mdsim::potentials::gpu::neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    boost::shared_ptr<box_type> box;
    boost::shared_ptr<potential_type> potential;
    boost::shared_ptr<force_type> force;
    boost::shared_ptr<neighbour_type> neighbour;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<host_potential_type> host_potential;

    modified_lennard_jones();
    void test();
};

template <typename float_type>
void modified_lennard_jones<float_type>::test()
{
    using namespace halmd::mdsim::gpu; // particle_kernel::{tagged, untagged}

    // place particles along the x-axis within one half of the box,
    // put every second particle at the origin
    unsigned int npart = particle->nbox;
    vector_type dx(0);
    dx[0] = box->edges()[0][0] / npart / 2;

    cuda::host::vector<float4> r_list(particle->g_r.size());
    for (unsigned int k = 0; k < r_list.size(); ++k) {
        vector_type r = (k % 2) ? k * dx : vector_type(0);
        unsigned int type = (k < particle->ntypes[0]) ? 0U : 1U;  // set particle type for a binary mixture
        r_list[k] = particle_kernel::tagged<vector_type>(r, type);
    }
    cuda::copy(r_list, particle->g_r);

    force->aux_enable();              // enable computation of auxiliary quantities
    force->compute();

    // read forces and other stuff from device
    cuda::host::vector<typename particle_type::gpu_vector_type> f_list(particle->g_f.size());
    cuda::copy(particle->g_f, f_list);

    cuda::host::vector<float> en_pot(force->potential_energy().size());
    cuda::copy(force->potential_energy(), en_pot);

    cuda::host::vector<float> hypervirial(force->hypervirial().size());
    cuda::copy(force->hypervirial(), hypervirial);

    const float_type tolerance = 20 * numeric_limits<float_type>::epsilon(); // FIXME the prefactor is an unjustified guess

    for (unsigned int i = 0; i < npart; ++i) {
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

        // FIXME the tolerance needs to cover both very large and vanishing forces
        BOOST_CHECK_SMALL(norm_inf(fval * r - f), max(norm_inf(fval * r), float_type(1)) * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 4 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, hypervirial[i], tolerance);
    }
}

template <typename float_type>
modified_lennard_jones<float_type>::modified_lennard_jones()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    vector<unsigned int> npart_list = list_of(1000)(2);
    float box_length = 100;
    float cutoff = box_length / 2;

    boost::array<float, 3> cutoff_array = list_of(cutoff)(cutoff)(cutoff);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(4.f);
    boost::array<unsigned, 3> index_m_array = list_of(12)(12)(24);
    boost::array<unsigned, 3> index_n_array = list_of(4)(2)(12);

    // create modules
    particle = boost::make_shared<particle_type>(npart_list);
    box = boost::make_shared<box_type>(particle->nbox, fixed_vector<double, dimension>(box_length));
    potential = boost::make_shared<potential_type>(
        particle->ntype, cutoff_array
      , epsilon_array, sigma_array, index_m_array, index_n_array
    );
    host_potential = boost::make_shared<host_potential_type>(
        particle->ntype, cutoff_array
      , epsilon_array, sigma_array, index_m_array, index_n_array
    );
    neighbour = boost::make_shared<neighbour_type>(particle);
    force = boost::make_shared<force_type>(potential, particle, box, neighbour);
}

BOOST_FIXTURE_TEST_CASE( modified_lennard_jones_gpu, device ) {
    modified_lennard_jones<float>().test();
}
#endif // WITH_CUDA
