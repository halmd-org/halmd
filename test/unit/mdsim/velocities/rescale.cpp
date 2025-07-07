/*
 * Copyright © 2025 Giorgia Marcelli
 * Copyright © 2011-2025 Felix Höfling
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

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/host/velocities/rescale.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/mdsim/gpu/velocities/rescale.hpp>
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
 * obtain particle data and verify that total energy equal target value
 */
template <typename particle_type>
void check_target_energies(particle_type const& particle, double target_energy, double tolerance)
{
    // obtain particle data: velocities, masses, potential energies
    std::vector<typename particle_type::velocity_type> velocity(particle.nparticle());
    BOOST_CHECK(
        get_velocity(particle, velocity.begin()) == velocity.end()
    );

    std::vector<typename particle_type::mass_type> mass(particle.nparticle());
    BOOST_CHECK(
        get_mass(particle, mass.begin()) == mass.end()
    );

    std::vector<typename particle_type::en_pot_type> potential_energy(particle.nparticle());
    BOOST_CHECK(
        particle.template get_data<typename particle_type::en_pot_type>(
            "potential_energy", potential_energy.begin()
        ) == potential_energy.end()
    );

    // check total energy for each particle
    for (unsigned int i = 0; i < particle.nparticle(); ++i) {
        auto const& v = velocity[i];
        double en_pot = potential_energy[i];
        double en_tot = en_pot + mass[i] * inner_prod(v, v) / 2;
//        BOOST_TEST_MESSAGE("v[" << i << "] = " << v << ", en_pot = " << en_pot);
        BOOST_CHECK_CLOSE_FRACTION(en_tot, target_energy, tolerance);
    }
}

/**
 * test energy rescaling using the rescale velocity module
 **/

template <typename modules_type>
struct rescale_test
{
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::boltzmann_type boltzmann_type;
    typedef typename modules_type::rescale_type rescale_type;
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

    void test();
    rescale_test();
};

template <typename modules_type>
void rescale_test<modules_type>::test()
{
    particle_group_type group(particle);

    // assign linear mass profile to particle array
    float_type scale_mass = float_type(0.2) / npart;
    set_mass(*particle, boost::make_transform_iterator(
        boost::make_counting_iterator(0)
      , [=](int i) { return float_type(0.9) + scale_mass * i; }
    ));

    // assign potential energy profile to particle array, from 0 to target_energy
    typedef typename particle_type::en_pot_type en_pot_type;
    float_type scale_energy = target_energy / (npart - 1);
    particle->template set_data<en_pot_type>(
        "potential_energy"
      , boost::make_transform_iterator(
            boost::make_counting_iterator(0)
          , [=](int i) {
                // ensure that values do not exceed the target value, be aware of rounding errors
                return std::min(scale_energy * i, static_cast<en_pot_type>(target_energy));
            }
    ));

    // assign initial velocities from Maxwell-Boltzmann distribution
    BOOST_TEST_MESSAGE("generate Maxwell-Boltzmann distribution");
    boltzmann->set();

    // compute initial energies
    double en_kin = get_mean_en_kin(*particle, group);
    double en_pot = get_mean_en_pot(*particle, group);
    double en_tot = en_kin + en_pot;

    BOOST_TEST_MESSAGE("energy before rescaling: total = " << en_tot
        << ", kinetic = " << en_kin << ", potential = " << en_pot);

    // apply velocity rescaling
    rescale_type rescaler(particle, target_energy, rescale_type::nve);
    rescaler.set();

    // test energy of each particle separately
    check_target_energies(*particle, target_energy, 10 * tolerance::value);

    // compute mean energies after rescaling
    en_kin = get_mean_en_kin(*particle, group);
    en_pot = get_mean_en_pot(*particle, group);
    en_tot = en_kin + en_pot;

    BOOST_TEST_MESSAGE("energy after rescaling: total = " << en_tot
        << ", kinetic = " << en_kin << ", potential = " << en_pot);

    // check that total energy is close to target value
    BOOST_CHECK_CLOSE_FRACTION(en_tot, target_energy, tolerance::value);

    // lower target energy and trap exception
    rescaler.set_target_energy(target_energy / 2);
    BOOST_CHECK_THROW(rescaler.set(), std::runtime_error);
}

template <typename modules_type>
rescale_test<modules_type>::rescale_test()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    npart = gpu ? 3000 : 30;
    temp = 2.2;
    target_energy = float(1.3);     // avoid possible round-off issues with lower precision

    // construct test modules, keep their dependencies only locally
    particle = std::make_shared<particle_type>(npart, 1);

    auto random = std::make_shared<random_type>();
    boltzmann = std::make_shared<boltzmann_type>(particle, random, temp);
}

// tolerance helper
template<typename float_type>
struct host_tolerance
{
    static constexpr double value = 2 * std::numeric_limits<float_type>::epsilon();
};

template<typename float_type>
constexpr double host_tolerance<float_type>::value;

// module traits
template <int dimension, typename float_type>
struct host_modules
{
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::host::particle_groups::all<particle_type> particle_group_type;

    typedef halmd::random::host::random random_type;
    typedef halmd::mdsim::host::velocities::boltzmann<dimension, float_type> boltzmann_type;
    typedef halmd::mdsim::host::velocities::rescale<dimension, float_type> rescale_type;

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

// dsfloat has effectively 43 bits, single float merely 24.
// For an unknown reason, rescale<dsfloat> seems to be not much more precise than using single precision,
// so let's mulitply by a huge, empirical factor
template<>
double const gpu_tolerance<halmd::dsfloat>::value = (1U << 17) * std::numeric_limits<float>::epsilon() / (1U << (43 - 24));

template<>
double const gpu_tolerance<float>::value = std::numeric_limits<float>::epsilon() / 2;

template <int dimension, typename float_type>
struct gpu_modules
{
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> particle_group_type;

    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef halmd::mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> boltzmann_type;
    typedef halmd::mdsim::gpu::velocities::rescale<dimension, float_type> rescale_type;

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
