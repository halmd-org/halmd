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
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/mdsim/gpu/velocity.hpp>
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

    std::shared_ptr<box_type> box;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<random_type> random;
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

    // compute initial energy
    double en_kin_before = thermo->en_kin();
    double en_pot_before = thermo->en_pot();
    double energy_before = en_kin_before + en_pot_before;

    BOOST_TEST_MESSAGE("Energy before rescale: total = " << energy_before
        << ", kinetic = " << en_kin_before
        << ", potential = " << en_pot_before);


    // apply rescale module
    rescale_type rescaler(particle, thermo, target_energy);
    rescaler.set();

    // Re-evaluate energies
    double en_kin_after = thermo->en_kin();
    double en_pot_after = thermo->en_pot();
    double energy_after = en_kin_after + en_pot_after;

    BOOST_TEST_MESSAGE("Energy after rescale: total = " << energy_after
        << ", kinetic = " << en_kin_after
        << ", potential = " << en_pot_after);

    // Check total energy is close to target
    float_type tolerance = 2 * std::numeric_limits<float_type>::epsilon();
    BOOST_TEST_MESSAGE("Target energy: " << target_energy);
    BOOST_CHECK_CLOSE_FRACTION(energy_after, target_energy, tolerance); // tolerance::value

    // Check center-of-mass velocity is small
    auto v_cm = get_v_cm(*particle, group);
    BOOST_TEST_MESSAGE("Center-of-mass velocity: " << v_cm);
    BOOST_CHECK_SMALL(norm_inf(v_cm), tolerance);
}

template <typename modules_type>
rescale_test<modules_type>::rescale_test()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    npart = gpu ? 10000 : 300;
    temp = 2.0;
    density = 0.3;
    target_energy = 1.0;

    double box_length = std::pow(npart / density, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i)
        edges(i, i) = box_length;

    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    random = std::make_shared<random_type>();
    boltzmann = std::make_shared<boltzmann_type>(particle, random, temp);
    // = std::make_shared<thermo_type>(particle);
    auto group = std::make_shared<particle_group_type>(particle);
    auto logger = std::make_shared<halmd::logger>();

    thermo = std::make_shared<thermo_type>(
        particle,
        group,
        box,
        [=]() { return std::pow(box_length, dimension); },
        logger
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

