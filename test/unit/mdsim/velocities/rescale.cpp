/*
 * Copyright © 2025 Giorgia Marcelli
 * Copyright © 2011-2017 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

// define boost unit test module and include set up 
#define BOOST_TEST_MODULE rescale
#include <boost/test/unit_test.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/host/velocities/rescale.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/mdsim/gpu/velocities/rescale.hpp>
# include <halmd/observables/gpu/thermodynamics.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/tools/cuda.hpp>
#endif
#include <test/tools/ctest.hpp>

#include <boost/assign.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <limits>

/**
 * test energy rescaling using the rescale velocity module
 **/

template <typename modules_type>
struct rescale_test
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::boltzmann_type boltzmann_type;
    typedef typename modules_type::rescale_type rescale_type;
    typedef typename modules_type::thermo_type thermo_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;
    static bool const gpu = modules_type::gpu;

    unsigned npart;
    double temp;
    double density;
    double target_energy;

    typedef typename modules_type::tolerance tolerance;

    std::shared_ptr<particle_type> particle;
    std::shared_ptr<boltzmann_type> boltzmann;
    std::shared_ptr<thermo_type> thermo;

    void test();
    rescale_test();
};

template <typename modules_type>
void rescale_test<modules_type>::test()
{
    particle_group_type group(particle);

    // assign mass profile
    float_type scale_mass = float_type(0.2) / npart;
    set_mass(*particle, boost::make_transform_iterator(
        boost::make_counting_iterator(0)
      , [=](int i) {
            return float_type(0.9) + scale_mass * i;
        }
    ));

    // assign initial velocities from Maxwell-Boltzmann distribution
    BOOST_TEST_MESSAGE("generate Maxwell-Boltzmann distribution");
    boltzmann->set();

    // output energies before rescaling
    BOOST_TEST_MESSAGE("Energy before rescaling: total = " << thermo->en_tot()
        << ", kinetic = " << thermo->en_kin()
        << ", potential = " << thermo->en_pot());

    // apply rescale module
    rescale_type rescaler(particle, target_energy);
    rescaler.set();

    // output energies after rescaling
    double energy_after = thermo->en_tot();
    BOOST_TEST_MESSAGE("Energy after rescaling: total = " << energy_after
        << ", kinetic = " << thermo->en_kin()
        << ", potential = " << thermo->en_pot());

    // check that total energy is close to target value
    float_type tolerance = 2 * std::numeric_limits<float_type>::epsilon();
    BOOST_CHECK_CLOSE_FRACTION(energy_after, target_energy, tolerance); // tolerance::value
}

template <typename modules_type>
rescale_test<modules_type>::rescale_test()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    npart = gpu ? 3000 : 30;
    temp = 2.2;
    density = 0.3;
    target_energy = 1.3;

    double volume = npart / density;
    double box_length = std::pow(volume, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i)
        edges(i, i) = box_length;

    // construct test modules, keep their dependencies only locally
    particle = std::make_shared<particle_type>(npart, 1);

    auto random = std::make_shared<random_type>();
    boltzmann = std::make_shared<boltzmann_type>(particle, random, temp);

    auto group = std::make_shared<particle_group_type>(particle);
    auto box = std::make_shared<box_type>(edges);
    auto logger = std::make_shared<halmd::logger>();
    thermo = std::make_shared<thermo_type>(
        particle, group, box
      , [=]() { return volume; }
      , logger
    );
}

// tolerance helper
template<typename float_type>
struct host_tolerance
{
    static constexpr double value = 2 * std::numeric_limits<float_type>::epsilon();
};

template struct host_tolerance<double>;

// module traits
template <int dimension, typename float_type>
struct host_modules
{
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::host::particle_groups::all<particle_type> particle_group_type;

    typedef halmd::random::host::random random_type;
    typedef halmd::mdsim::host::velocities::boltzmann<dimension, float_type> boltzmann_type;
    typedef halmd::mdsim::host::velocities::rescale<dimension, float_type> rescale_type;
    typedef halmd::observables::host::thermodynamics<dimension, float_type> thermo_type;

    static bool const gpu = false;
    typedef host_tolerance<float_type> tolerance;
};

// Register tests with Boost.Test
#ifndef USE_HOST_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE(rescale_host_2d)
{
    rescale_test<host_modules<2, double>>().test();
}
BOOST_AUTO_TEST_CASE(rescale_host_3d)
{
    rescale_test<host_modules<3, double>>().test();
}
#else
BOOST_AUTO_TEST_CASE(rescale_host_2d)
{
    rescale_test<host_modules<2, float>>().test();
}
BOOST_AUTO_TEST_CASE(rescale_host_3d)
{
    rescale_test<host_modules<3, float>>().test();
}
#endif

#ifdef HALMD_WITH_GPU
template<typename T>
struct gpu_tolerance
{
    static double const value;
};

// dsfloat has effectively 43 bits, single float merely 24,
// multiply by number of particles to get a sharp upper bound
template<>
double const gpu_tolerance<halmd::dsfloat>::value = 10000 * std::numeric_limits<float>::epsilon() / (1U << (43 - 24));

template<>
double const gpu_tolerance<float>::value = std::numeric_limits<float>::epsilon() / 4;    // yields 0.5 ulp, the 4 is empirical

template <int dimension, typename float_type>
struct gpu_modules
{
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> particle_group_type;

    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef halmd::mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> boltzmann_type;
    typedef halmd::mdsim::gpu::velocities::rescale<dimension, float_type> rescale_type;
    typedef halmd::observables::gpu::thermodynamics<dimension, float_type> thermo_type;

    static bool const gpu = true;
    typedef gpu_tolerance<float_type> tolerance;
};

# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( rescale_gpu_float_2d, set_cuda_device ) {
    rescale_test<gpu_modules<2, float>>().test();
}
BOOST_FIXTURE_TEST_CASE( rescale_gpu_float_3d, set_cuda_device ) {
    rescale_test<gpu_modules<3, float>>().test();
}
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( rescale_gpu_dsfloat_2d, set_cuda_device ) {
    rescale_test<gpu_modules<2, halmd::dsfloat>>().test();
}
BOOST_FIXTURE_TEST_CASE( rescale_gpu_dsfloat_3d, set_cuda_device ) {
    rescale_test<gpu_modules<3, halmd::dsfloat>>().test();
}
# endif
#endif // HALMD_WITH_GPU
