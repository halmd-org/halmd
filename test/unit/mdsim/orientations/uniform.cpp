/*
 * Copyright © 2023 Jaslo Ziska
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

#define BOOST_TEST_MODULE uniform
#include <boost/test/unit_test.hpp>

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/orientations/uniform.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/orientations/uniform.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

/**
 * test initialisation of particle orientations: uniform module
 */

template <typename modules_type>
struct uniform
{
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::orientation_type orientation_type;
    typedef typename particle_type::vector_type::value_type float_type;
    static unsigned int const dimension = particle_type::vector_type::static_size;
    static bool const gpu = modules_type::gpu;

    unsigned npart;

    std::shared_ptr<particle_type> particle;
    std::shared_ptr<random_type> random;
    std::shared_ptr<orientation_type> orientation;

    void test();
    uniform();
};

template <typename modules_type>
void uniform<modules_type>::test()
{
    // generate orientation distribution
    BOOST_TEST_MESSAGE("generate uniform orientation distribution");
    orientation->set();

    // typename particle_type::orientation_array_type const& orientations = read_cache(particle->orientation());
    std::vector<typename particle_type::orientation_type> orientations(npart);
    BOOST_CHECK(get_orientation(*particle, orientations.begin()) == orientations.end());

    halmd::accumulator<typename particle_type::orientation_type> a;
    for (unsigned int i = 0; i < npart; ++i) {
        BOOST_CHECK_CLOSE(norm_2(orientations[i]), 1, 1e-3); // FIXME: what tolerance to use?
        a(orientations[i]);
    }

    // check the mean of the individual orientations
    // mean = 0, variance = 1/dimension
    BOOST_CHECK_EQUAL(count(a), npart);
    float_type tol = 4.5 / std::sqrt(dimension * (npart - 1.));
    BOOST_CHECK_SMALL(mean(a)[0], tol);
    BOOST_CHECK_SMALL(mean(a)[1], tol);
    if (dimension == 3) BOOST_CHECK_SMALL(mean(a)[2], tol);

    // check the distance of the "random walk" of all orientations
    // mean = 0, variance = npart
    tol = 4.5 * std::sqrt(npart);
    BOOST_CHECK_SMALL(norm_2(sum(a)), tol);
}

template <typename modules_type>
uniform<modules_type>::uniform()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    npart = gpu ? 10000 : 3000;
    particle = std::make_shared<particle_type>(npart, 1);
    random = std::make_shared<random_type>();
    orientation = std::make_shared<orientation_type>(particle, random);
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::random::host::random random_type;
    typedef halmd::mdsim::host::orientations::uniform<dimension, float_type> orientation_type;
    static bool const gpu = false;
};

#ifndef USE_HOST_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE(uniform_host_2d) {
    uniform<host_modules<2, double>>().test();
}
BOOST_AUTO_TEST_CASE(uniform_host_3d) {
    uniform<host_modules<3, double>>().test();
}
#else
BOOST_AUTO_TEST_CASE(uniform_host_2d) {
    uniform<host_modules<2, float>>().test();
}
BOOST_AUTO_TEST_CASE(uniform_host_3d) {
    uniform<host_modules<3, float>>().test();
}
#endif

#ifdef HALMD_WITH_GPU
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef halmd::mdsim::gpu::orientations::uniform<dimension, float_type, halmd::random::gpu::rand48> orientation_type;
    static bool const gpu = true;
};

# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE(uniform_gpu_float_2d, halmd::device) {
    uniform<gpu_modules<2, float>>().test();
}
BOOST_FIXTURE_TEST_CASE(uniform_gpu_float_3d, halmd::device) {
    uniform<gpu_modules<3, float>>().test();
}
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE(uniform_gpu_dsfloat_2d, halmd::device) {
    uniform<gpu_modules<2, halmd::dsfloat>>().test();
}
BOOST_FIXTURE_TEST_CASE(uniform_gpu_dsfloat_3d, halmd::device) {
    uniform<gpu_modules<3, halmd::dsfloat>>().test();
}
# endif
#endif // HALMD_WITH_GPU
