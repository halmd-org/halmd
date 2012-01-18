/*
 * Copyright © 2011  Michael Kopp
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

/** test power law Potential with divergence at a finite distance (core)
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
    unsigned int ntype = 2;  // test a binary mixture
    boost::array<float, 3> cutoff_array = list_of(5.f)(6.f)(7.f);
    boost::array<float, 3> core_array = list_of(0.4f)(0.5f)(0.6f);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(2.2f);
    boost::array<unsigned, 3> index_array = list_of(6u)(12u)(7u);

    // construct module
    potential_type potential(ntype, cutoff_array, core_array, epsilon_array, sigma_array, index_array);

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

    matrix_type cutoff = potential.r_cut_sigma();
    BOOST_CHECK(cutoff(0, 0) == cutoff_array[0]);
    BOOST_CHECK(cutoff(0, 1) == cutoff_array[1]);
    BOOST_CHECK(cutoff(1, 0) == cutoff_array[1]);
    BOOST_CHECK(cutoff(1, 1) == cutoff_array[2]);

    matrix_type core = potential.r_core_sigma();
    BOOST_CHECK(core(0, 0) == core_array[0]);
    BOOST_CHECK(core(0, 1) == core_array[1]);
    BOOST_CHECK(core(1, 0) == core_array[1]);
    BOOST_CHECK(core(1, 1) == core_array[2]);

    uint_matrix_type index = potential.index();
    BOOST_CHECK(index(0, 0) == index_array[0]);
    BOOST_CHECK(index(0, 1) == index_array[1]);
    BOOST_CHECK(index(1, 0) == index_array[1]);
    BOOST_CHECK(index(1, 1) == index_array[2]);

    // evaluate some points of the potential, force, and hypervirial
    typedef boost::tuple<double, double, double, double> tuple_type;
    const double eps = numeric_limits<double>::epsilon();
    double tolerance = 1; // dummy

    // expected results (r, fval, en_pot, hvir)
    // interaction AA: ε=1, σ=1, rc=5σ, r_core=0.4σ, n=6
    vector<tuple_type> results = tuple_list_of
        (5.00000000000000e-01,  1.19999999987334e+08,  9.99999999894453e+05,  1.01999999989234e+09)
        (1.00000000000000e+00,  2.14333649588151e+02,  2.14333649588151e+01,  2.28622559560695e+03)
        (3.00000000000000e+00,  2.40890735993852e-03,  3.13157956792008e-03,  1.53428868771469e-01)
        (5.00000000000000e+00, -7.07088373360112e-21, -2.71050543121376e-20, -1.16823296468192e-18)
        (1.00000000000000e+01, -6.51694969216805e-06, -1.04271195074689e-04, -4.10024751465573e-03);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 0);  // interaction AA

        // tolerance for temp = (r-r_core)^-(n+2)
        double tolerance_ = eps * index_array[0] * ( 1 + max(get<0>(a), double(cutoff_array[0]))/(get<0>(a) - cutoff_array[0]) );
        // potential without cutoff-correction
        double pot_nocut = epsilon_array[0] * pow( get<0>(a)/sigma_array[0] - core_array[0] , -index_array[0] );
        // cutoff energy
        double V0 = epsilon_array[0] * pow( cutoff_array[0] - core_array[0] , -index_array[0] );

        BOOST_CHECK_CLOSE(en_pot, get<2>(a), tolerance_ + eps * (1 + max(pot_nocut, V0)/abs(pot_nocut - V0))); // errors due to subtraction and multiplication
        double tolerance_fval = tolerance_ + eps * (4 + max(get<0>(a), double(cutoff_array[0]))/abs(get<0>(a)-cutoff_array[0]));
        BOOST_CHECK_CLOSE(fval, get<1>(a), tolerance_fval);
        BOOST_CHECK_CLOSE(hvir, get<3>(a), tolerance_fval + eps * 5);
    };

    // interaction AB: ε=0.5, σ=2, rc=6σ, r_core=0.5σ, n=12
    results = tuple_list_of
        (1.50000000000000e+00,  1.34217728000000e+08,  8.38860799999999e+06,  1.14756157440000e+10)
        (3.00000000000000e+00,  9.99999985497268e-01,  4.99999992748634e-01,  1.66499997585295e+02)
        (5.00000000000000e+00,  7.32378366802688e-05,  1.22063061133781e-04,  2.79219252343525e-02)
        (9.30000000000000e+00,  1.85130452966199e-09,  1.19085163870507e-08,  2.17222125536800e-06)
        (1.50000000000000e+01, -4.12299566330883e-10, -7.21524241079045e-09, -1.19934998848752e-06);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 0, 1);  // interaction AB
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };

    // interaction BB: ε=0.25, σ=2.2, rc=7σ, r_core=0.6σ, n=7
    results = tuple_list_of
        (1.75000000000000e+00,  2.13408208954871e+05,  2.29413824626486e+04,  2.06252209855177e+07)
        (4.00000000000000e+00,  4.10019316439325e-02,  6.27915296032796e-02,  7.17717394507524e+00)
        (7.00000000000000e+00,  5.61802960847151e-05,  3.19104081761182e-04,  2.43877874031972e-02)
        (1.22000000000000e+01, -2.30709843967200e-07, -4.37478597840433e-06, -2.73700859138740e-04)
        (1.90000000000000e+01, -1.60760926012554e-07, -7.71468718087676e-06, -4.40906116352052e-04);

    BOOST_FOREACH (tuple_type a, results) {
        double rr = std::pow(get<0>(a), 2);
        double fval, en_pot, hvir;
        tie(fval, en_pot, hvir) = potential(rr, 1, 1);  // interaction BB
        BOOST_CHECK_CLOSE_FRACTION(fval, get<1>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(en_pot, get<2>(a), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(hvir, get<3>(a), tolerance);
    };
}

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

    power_law_with_core();
    void test();
};

template <typename float_type>
void power_law_with_core<float_type>::test()
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
//        r_list[k] = particle_kernel::tagged<vector_type>(r, type); //< GCC 4.4.1 crashes with internal compiler error
        r_list[k] = r;
        union { float f; uint32_t i; } u; u.i = type; r_list[k].w = u.f;
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

    const float_type tolerance = 10 * numeric_limits<float_type>::epsilon();

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
power_law_with_core<float_type>::power_law_with_core()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    vector<unsigned int> npart_list = list_of(1000)(2);
    float box_length = 100;
    float cutoff = box_length / 2;

    boost::array<float, 3> cutoff_array = list_of(cutoff)(cutoff)(cutoff);
    boost::array<float, 3> core_array = list_of(0.4f)(0.5f)(0.6f);
    boost::array<float, 3> epsilon_array = list_of(1.f)(.5f)(.25f);
    boost::array<float, 3> sigma_array = list_of(1.f)(2.f)(2.2f);
    boost::array<unsigned, 3> index_array = list_of(6u)(12u)(7u);

    // create modules
    particle = make_shared<particle_type>(npart_list);
    box = make_shared<box_type>(particle->nbox, fixed_vector<double, dimension>(box_length));
    potential = make_shared<potential_type>(particle->ntype, cutoff_array, core_array, epsilon_array, sigma_array, index_array);
    host_potential = make_shared<host_potential_type>(particle->ntype, cutoff_array, core_array, epsilon_array, sigma_array,index_array);
    neighbour = make_shared<neighbour_type>(particle);
    force = make_shared<force_type>(potential, particle, box, neighbour);
}

#ifdef WITH_CUDA
BOOST_FIXTURE_TEST_CASE( power_law_with_core_gpu, device ) {
    power_law_with_core<float>().test();
}
#endif // WITH_CUDA
