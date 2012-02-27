/*
 * Copyright © 2011-2012  Felix Höfling
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

#define BOOST_TEST_MODULE power_law
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <cmath> // std::pow
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/power_law.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_kernel.cuh>
# include <halmd/mdsim/gpu/potentials/power_law.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/gpu/neighbour_chain.hpp>
#endif


using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace std;

/** test power law potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the pair_trunc force module in two dimensions to
 *  compute some values of the potential which are compared against the host
 *  module. This requires a special neighbour list module with only one defined
 *  neighbour per particle.
 */

BOOST_AUTO_TEST_CASE( power_law_host )
{
    typedef mdsim::host::potentials::power_law<double> potential_type;
    typedef potential_type::matrix_type matrix_type;
    typedef potential_type::uint_matrix_type uint_matrix_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    boost::array<float, 3> cutoff_array = list_of(5.f)(5.f)(5.f);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(4.f);
    boost::array<unsigned, 3> index_array = list_of(12)(24)(3);

    // construct module
    potential_type potential(ntype, cutoff_array, epsilon_array, sigma_array, index_array);

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

    uint_matrix_type index = potential.index();
    BOOST_CHECK(index(0, 0) == index_array[0]);
    BOOST_CHECK(index(0, 1) == index_array[1]);
    BOOST_CHECK(index(1, 0) == index_array[1]);
    BOOST_CHECK(index(1, 1) == index_array[2]);

    // evaluate some points of the potential, force, and hypervirial
    typedef boost::tuple<double, double, double, double> tuple_type;
    double const eps = numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot, hvir) for ε=1, σ=1, n=12, rc=5σ
    vector<tuple_type> results = tuple_list_of
        (0.2, 7.32421875e10, 2.44140625e8, 3.515625e10)
        (0.5, 196608., 4095.999999995904, 589824.)
        (1., 12., 0.999999995904, 144.)
        (2., 0.000732421875, 0.000244136529, 0.03515625)
        (10., 1.2e-13, -4.095e-9, 1.44e-10);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 0);  // interaction AA

        double tolerance = (1 + index_array[0] / 2) * eps;
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction AB: ε=.5, σ=2, n=24, rc=5σ
    results = tuple_list_of
        (0.2, 3.e26, 5.e23, 2.88e26)
        (0.5, 1.3510798882111488e16,  1.40737488355328e14, 8.106479329266893e16)
        (1., 2.01326592e8, 8.388608e6, 4.831838208e9)
        (2., 3., 0.5, 288.)
        (10., 2.01326592e-18, 0., 4.831838208e-15);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 1);  // interaction AB

        double tolerance = (1 + index_array[1] / 2) * eps;
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction BB: ε=.25, σ=4, n=3, rc=5σ
    results = tuple_list_of
        (0.2, 150000., 1999.998, 18000.)
        (0.5, 1536., 127.998, 1152.)
        (1., 48., 15.998, 144.)
        (2., 1.5, 1.998, 18.)
        (10., 0.00048, 0.014, 0.144);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 1, 1);  // interaction BB

        double tolerance = (1 + index_array[2] / 2) * eps;
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };
}

#ifdef WITH_CUDA

template <typename float_type>
struct power_law
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::power_law<float_type> potential_type;
    typedef mdsim::host::potentials::power_law<double> host_potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef test::unit::mdsim::potentials::gpu::neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    shared_ptr<box_type> box;
    shared_ptr<potential_type> potential;
    shared_ptr<force_type> force;
    shared_ptr<neighbour_type> neighbour;
    shared_ptr<particle_type> particle;
    shared_ptr<host_potential_type> host_potential;

    power_law();
    void test();
};

template <typename float_type>
void power_law<float_type>::test()
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

        // estimated upper bound on floating-point error for power function
        float_type const eps = numeric_limits<float_type>::epsilon();
        unsigned int index = host_potential->index()(type1, type2);
        float_type tolerance = index * eps;

        // check both absolute and relative error
        BOOST_CHECK_SMALL(norm_inf(fval * r - f), max(norm_inf(fval * r), float_type(1)) * tolerance);

        // the prefactor is not justified, it is needed for absolute differences on the order 1e-18
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 5 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, hypervirial[i], tolerance);
    }
}

template <typename float_type>
power_law<float_type>::power_law()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    vector<unsigned int> npart_list = list_of(1000)(2);
    float box_length = 100;
    float cutoff = box_length / 2;

    boost::array<float, 3> cutoff_array = list_of(cutoff)(cutoff)(cutoff);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(4.f);
    boost::array<unsigned, 3> index_array = list_of(6)(12)(24);

    // create modules
    particle = make_shared<particle_type>(npart_list);
    box = make_shared<box_type>(particle->nbox, fixed_vector<double, dimension>(box_length));
    potential = make_shared<potential_type>(
        particle->ntype, cutoff_array
      , epsilon_array, sigma_array, index_array
    );
    host_potential = make_shared<host_potential_type>(
        particle->ntype, cutoff_array
      , epsilon_array, sigma_array, index_array
    );
    neighbour = make_shared<neighbour_type>(particle);
    force = make_shared<force_type>(potential, particle, box, neighbour);
}

BOOST_FIXTURE_TEST_CASE( power_law_gpu, device ) {
    power_law<float>().test();
}
#endif // WITH_CUDA
