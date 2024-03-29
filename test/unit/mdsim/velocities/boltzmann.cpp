/*
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

#define BOOST_TEST_MODULE boltzmann
#include <boost/test/unit_test.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/random/host/random.hpp>
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
 * test initialisation of particle velocities: boltzmann module
 */

template <typename modules_type>
struct boltzmann
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::velocity_type velocity_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;
    static bool const gpu = modules_type::gpu;

    unsigned npart;
    double temp;
    double density;

    typedef typename modules_type::tolerance tolerance;

    std::shared_ptr<box_type> box;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<random_type> random;
    std::shared_ptr<velocity_type> velocity;

    void test();
    boltzmann();
};

template <typename modules_type>
void boltzmann<modules_type>::test()
{
    particle_group_type group(particle);

    // assign a flat distribution of particle masses
    float_type scale_mass = float_type(0.2) / npart;
    set_mass(*particle, boost::make_transform_iterator(
        boost::make_counting_iterator(0)
      , [=](int i) {
            return float_type(0.9) + scale_mass * i;
        }
    ));

    // generate velocity distribution
    BOOST_TEST_MESSAGE("generate Maxwell-Boltzmann distribution");
    velocity->set();

    //
    // test velocity distribution of final state
    //
    halmd::fixed_vector<double, dimension> v_cm = get_v_cm(*particle, group);
    double en_kin = get_mean_en_kin(*particle, group);
    BOOST_TEST_MESSAGE("Centre-of-mass velocity: " << npart << " " << v_cm);
    BOOST_TEST_MESSAGE("Instantaneous temperature: " << 2 * en_kin / dimension);

    // centre-of-mass velocity ⇒ mean of velocity distribution
    // each particle is an independent "measurement",
    // tolerance is 4.5σ, σ = √(<v_x²> / (N - 1)) where <v_x²> = k T / m,
    // with this choice, a single test passes with 99.999% probability
    double vcm_tolerance = 4.5 * sqrt(temp / (npart - 1));
    BOOST_TEST_MESSAGE("Absolute tolerance on instantaneous centre-of-mass velocity: " << vcm_tolerance);
    BOOST_CHECK_SMALL(norm_inf(get_v_cm(*particle, group)), vcm_tolerance);  //< norm_inf tests the max. value

    // temperature ⇒ use distribution of molecular kinetic energies from Maxwell—Boltzmann
    // distribution, each degree of freedom is an independent "measurement",
    // tolerance is 4.5σ, σ = √<ΔT²> where <ΔT²> / T² = 2 / (dimension × N - 1)
    double rel_temp_tolerance = 4.5 * sqrt(2. / (dimension * npart - 1));
    BOOST_TEST_MESSAGE("Relative tolerance on instantaneous temperature: " << rel_temp_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(2 * en_kin / dimension, temp, rel_temp_tolerance);

    //
    // test shifting and rescaling
    //
    // multiplication of the velocities by a constant factor
    double scale = 1.5;
    rescale_velocity_group(*particle, group, scale);
    BOOST_CHECK_CLOSE_FRACTION(2 * get_mean_en_kin(*particle, group) / dimension, scale * scale * temp, rel_temp_tolerance);

    // shift mean velocity to zero
    v_cm = get_v_cm(*particle, group);
    shift_velocity_group(*particle, group, -v_cm);
    vcm_tolerance = modules_type::tolerance::value;
    BOOST_CHECK_SMALL(norm_inf(get_v_cm(*particle, group)), vcm_tolerance);

    // first shift, then rescale in one step
    shift_rescale_velocity_group(*particle, group, v_cm, 1 / scale);
    BOOST_CHECK_CLOSE_FRACTION(2 * get_mean_en_kin(*particle, group) / dimension, temp, rel_temp_tolerance);
    BOOST_CHECK_SMALL(norm_inf(get_v_cm(*particle, group) - v_cm), vcm_tolerance);
}

template <typename modules_type>
boltzmann<modules_type>::boltzmann()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    npart = gpu ? 10000 : 3000;
    temp = 2.0;
    density = 0.3;
    double box_length = std::pow(npart / density, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }

    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    random = std::make_shared<random_type>();
    velocity = std::make_shared<velocity_type>(particle, random, temp);
}

template<typename float_type>
struct host_tolerance
{
    static constexpr double value = 2 * std::numeric_limits<float_type>::epsilon();
};

template <int dimension, typename float_type>
struct host_modules
{
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::random::host::random random_type;
    typedef halmd::mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    static bool const gpu = false;
    typedef host_tolerance<float_type> tolerance;
};

#ifndef USE_HOST_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE( boltzmann_host_2d ) {
    boltzmann<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( boltzmann_host_3d ) {
    boltzmann<host_modules<3, double> >().test();
}
#else
BOOST_AUTO_TEST_CASE( boltzmann_host_2d ) {
    boltzmann<host_modules<2, float> >().test();
}
BOOST_AUTO_TEST_CASE( boltzmann_host_3d ) {
    boltzmann<host_modules<3, float> >().test();
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
    typedef halmd::mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    static bool const gpu = true;
    typedef gpu_tolerance<float_type> tolerance;
};

# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( boltzmann_gpu_float_2d, set_cuda_device ) {
    boltzmann<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( boltzmann_gpu_float_3d, set_cuda_device ) {
    boltzmann<gpu_modules<3, float> >().test();
}
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( boltzmann_gpu_dsfloat_2d, set_cuda_device ) {
    boltzmann<gpu_modules<2, halmd::dsfloat> >().test();
}
BOOST_FIXTURE_TEST_CASE( boltzmann_gpu_dsfloat_3d, set_cuda_device ) {
    boltzmann<gpu_modules<3, halmd::dsfloat> >().test();
}
# endif
#endif // HALMD_WITH_GPU
