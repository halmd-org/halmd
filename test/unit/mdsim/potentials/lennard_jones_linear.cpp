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

#define BOOST_TEST_MODULE lennard_jones_linear
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <cmath> // std::pow
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/lennard_jones_linear.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_kernel.cuh>
# include <halmd/mdsim/gpu/potentials/lennard_jones_linear.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/mdsim/potentials/gpu/neighbour_chain.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace boost::assign; // list_of
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

BOOST_AUTO_TEST_CASE( lennard_jones_linear_host )
{
    typedef mdsim::host::potentials::lennard_jones_linear<double> potential_type;
    typedef potential_type::matrix_type matrix_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    boost::array<float, 3> cutoff_array = list_of(5.f)(5.f)(5.f);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(4.f);

    // construct module
    potential_type potential(ntype, cutoff_array, epsilon_array, sigma_array);

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

    // evaluate some points of the potential, force, and hypervirial
    typedef boost::tuple<double, double, double, double> tuple_type;
    const double tolerance = 5 * numeric_limits<double>::epsilon();

    // expected results (r, fval, en_pot, hvir) for ε=1, σ=1, rc=5σ
    vector<tuple_type> results = tuple_list_of
        (0.2, 2.929593750000015e11, 9.765000000017304e8, 1.406227499999999e11)
        (0.5, 780288.0006143214, 16128.00163820667, 2.35007999984642e6)
        (1., 24.0003071606784, 0.0014846263296, 431.9996928393216)
        (2., -0.0906667321608, -0.0603459718488, -2.1099893213568)
        (10., 0.00003047606832, -0.001283819772, -0.003215606208);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 0);  // interaction AA
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction AB: ε=.5, σ=2, rc=5σ
    results = tuple_list_of
        (0.2, 5.999997e14, 1.999998000000001e12, 2.87999928e14)
        (0.5, 1.610416128000154e9, 3.35462400008575e7, 4.831543295999962e9)
        (1., 97536.00007679017, 8064.000819103334, 1.17503999992321e6)
        (2., 3.0000383950848, 0.0007423131648, 215.9998464196608)
        (10., 0, 0, -0.005374722048);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 1);  // interaction AB
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction BB: ε=.25, σ=4, rc=5σ
    results = tuple_list_of
        (0.2, 1.2287999904e18, 4.095999936e15, 5.89823997696e17)
        (0.5, 3.298528591872e12, 6.871921459200044e10, 9.8955952128e12)
        (1., 2.013020160000192e8, 1.677312000042875e7, 2.415771647999981e9)
        (2., 12192.00000959877, 4032.000409551667, 587519.9999616049)
        (10., -0.00024182697984, -0.003823251456, -0.14523205632);

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
struct lennard_jones_linear
{
    enum { dimension = 2 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::lennard_jones_linear<float_type> potential_type;
    typedef mdsim::host::potentials::lennard_jones_linear<double> host_potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef test::unit::mdsim::potentials::gpu::neighbour_chain<dimension, float_type> neighbour_type;

    typedef typename particle_type::vector_type vector_type;

    shared_ptr<box_type> box;
    shared_ptr<potential_type> potential;
    shared_ptr<force_type> force;
    shared_ptr<neighbour_type> neighbour;
    shared_ptr<particle_type> particle;
    shared_ptr<host_potential_type> host_potential;

    lennard_jones_linear();
    void test();
};

template <typename float_type>
void lennard_jones_linear<float_type>::test()
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

    // FIXME find better estimate for floating point error, a global estimate
    // seems unsuitable due to the subtractions in the potential computation
    const float_type tolerance = 5e2 * 10 * numeric_limits<float_type>::epsilon();

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

        BOOST_CHECK_SMALL(norm_inf(fval * r - f), norm_inf(f) * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot_, en_pot[i], 4 * tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, hypervirial[i], tolerance);
    }
}

template <typename float_type>
lennard_jones_linear<float_type>::lennard_jones_linear()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    vector<unsigned int> npart_list = list_of(1000)(2);
    float box_length = 100;
    float cutoff = box_length / 2;

    boost::array<float, 3> cutoff_array = list_of(cutoff)(cutoff)(cutoff);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(4.f);

    // create modules
    particle = make_shared<particle_type>(npart_list);
    box = make_shared<box_type>(particle->nbox, fixed_vector<double, dimension>(box_length));
    potential = make_shared<potential_type>(particle->ntype, cutoff_array, epsilon_array, sigma_array);
    host_potential = make_shared<host_potential_type>(particle->ntype, cutoff_array, epsilon_array, sigma_array);
    neighbour = make_shared<neighbour_type>(particle);
    force = make_shared<force_type>(potential, particle, box, neighbour);
}

BOOST_FIXTURE_TEST_CASE( lennard_jones_linear_gpu, device ) {
    lennard_jones_linear<float>().test();
}
#endif // WITH_CUDA
