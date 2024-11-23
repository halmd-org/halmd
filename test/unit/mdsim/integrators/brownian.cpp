/*
 * Copyright © 2023 Jaslo Ziska
 * Copyright © 2015 Manuel Dibak
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

#define BOOST_TEST_MODULE brownian
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <functional>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/forces/external.hpp>
#include <halmd/mdsim/host/integrators/brownian.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/potentials/external/harmonic.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/dynamics/blocking_scheme.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/host/dynamics/mean_quartic_displacement.hpp>
#include <halmd/observables/host/dynamics/mean_square_displacement.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
# include <halmd/algorithm/gpu/apply_kernel.hpp>
# include <halmd/mdsim/gpu/forces/external.hpp>
# include <halmd/mdsim/gpu/integrators/brownian.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/potentials/external/harmonic.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/observables/gpu/dynamics/mean_quartic_displacement.hpp>
# include <halmd/observables/gpu/dynamics/mean_square_displacement.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace halmd;
using namespace std;

template <typename modules_type>
struct brownian_free
{
    typedef mdsim::clock clock_type;

    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::integrator_type integrator_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::phase_space_type phase_space_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::sample_type sample_type;

    typedef typename modules_type::msd_type msd_type;
    typedef typename modules_type::mqd_type mqd_type;

    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;

    static bool const gpu = modules_type::gpu;

    size_t steps;
    double density;
    double temperature;
    double timestep;
    typename integrator_type::matrix_type diff_const;
    double maximum_lag_time;
    double resolution;
    unsigned int block_size;
    unsigned int npart;
    fixed_vector<double, dimension> box_ratios;
    fixed_vector<double, dimension> slab;

    std::shared_ptr<box_type> box;
    std::shared_ptr<clock_type> clock;
    std::shared_ptr<integrator_type> integrator;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<particle_group_type> particle_group;
    std::shared_ptr<phase_space_type> phase_space;
    std::shared_ptr<position_type> position;
    std::shared_ptr<random_type> random;

    brownian_free();
    void test();
};

/** solve the differential equation @f$ \dot r = v = const. @f$ */
template <typename modules_type>
void brownian_free<modules_type>::test()
{
    // construct blocking scheme module (contains the logic)
    observables::dynamics::blocking_scheme blocking_scheme(
        clock, maximum_lag_time, resolution, block_size, 1, 100
    );

    // allocate space for block samples
    auto block_sample = std::make_shared<observables::samples::blocking_scheme<sample_type>>(
        [=]() { return phase_space->template acquire<sample_type>("position"); }
      , blocking_scheme.count()
      , blocking_scheme.block_size()
    );

    // construct correlation functions and connect to blocking scheme logic
    auto correlation_msd = make_shared<observables::dynamics::correlation<msd_type>>(
        make_shared<msd_type>(), block_sample, block_sample
    );
    blocking_scheme.on_correlate(correlation_msd);

    auto correlation_mqd = make_shared<observables::dynamics::correlation<mqd_type>>(
        make_shared<mqd_type>(), block_sample, block_sample
    );
    blocking_scheme.on_correlate(correlation_mqd);

    blocking_scheme.on_sample(block_sample);

    BOOST_TEST_MESSAGE("run Brownian integrator over " << steps << " steps");
    for (size_t i = 0; i < steps; ++i) {
        integrator->integrate();
        clock->advance();
        blocking_scheme.sample();
    }
    blocking_scheme.finalise();

    auto time = blocking_scheme.time()[0];
    auto msd = correlation_msd->result()[0];
    auto mqd = correlation_mqd->result()[0];

    for (size_t i = 0; i < size_t(maximum_lag_time / timestep); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(mean(msd[i]), 2 * dimension * diff_const(0, 0) * time[i], 4.5 * abs(error_of_mean(msd[i])));
        // BOOST_CHECK_CLOSE_FRACTION(mean(mqd[i]), 60 * time[i] * time[i], 4.5 * 120 * time[i] * error_of_mean(mqd[i])); // TODO: value for 2D
    }
}

/**
 * Initialize integrator and dependencies, set basic parameters.
 */
template <typename modules_type>
brownian_free<modules_type>::brownian_free()
{
    typedef fixed_vector<double, dimension> vector_type;

    timestep = 0.01;
    maximum_lag_time = 7;

    // run for as many steps as possible, wrap around the box for about 10 times
    // adjusted so total simulation length is constant when timestep changes
    steps = (gpu ? 1 : 100) * maximum_lag_time / timestep;
    resolution = 0.01;
    block_size = 10000;

    // a low density implies large values of the position vectors
    density = 1e-6;
    temperature = 1;
    // optimize filling of fcc lattice, use only few particles on the host
    npart = gpu ? 2500 : 5;

    box_ratios = (dimension == 3) ? vector_type{1., 2., 1.01} : vector_type{1., 2.};
    double det = accumulate(box_ratios.begin(), box_ratios.end(), 1., multiplies<double>());
    double volume = npart / density;
    double edge_length = pow(volume / det, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = edge_length * box_ratios[i];
    }
    slab = 1;

    //  DIFFUSION CONSTANT
    //  @param1 displacement
    //  @param2 rotational
    diff_const = typename integrator_type::matrix_type(1, 2);
    diff_const <<= 1.0, 1.0;

    particle = std::make_shared<particle_type>(npart, 1);
    particle_group = std::make_shared<particle_group_type>(particle);
    box = std::make_shared<box_type>(edges);
    random = std::make_shared<random_type>();
    integrator = std::make_shared<integrator_type>(
        particle, random, box, timestep, temperature, diff_const, integrator_type::integrate_both
    );
    position = std::make_shared<position_type>(particle, box, slab);
    clock = std::make_shared<clock_type>();
    clock->set_timestep(integrator->timestep());
    phase_space = std::make_shared<phase_space_type>(particle, particle_group, box);

    // set positions
    position->set();
}

template <typename modules_type>
struct brownian_harmonic
{
    typedef mdsim::clock clock_type;

    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::potential_type potential_type;
    typedef typename modules_type::force_type force_type;
    typedef typename modules_type::integrator_type integrator_type;
    typedef typename modules_type::phase_space_type phase_space_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::sample_type sample_type;

    typedef typename modules_type::msd_type msd_type;
    typedef typename modules_type::mqd_type mqd_type;

    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;

    static bool const gpu = modules_type::gpu;

    size_t steps;
    double temperature;
    double timestep;
    typename integrator_type::matrix_type diff_const;
    typename potential_type::scalar_container_type stiffness;
    double maximum_lag_time;
    double resolution;
    unsigned int block_size;
    unsigned int npart;

    std::shared_ptr<box_type> box;
    std::shared_ptr<clock_type> clock;
    std::shared_ptr<force_type> force;
    std::shared_ptr<integrator_type> integrator;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<particle_group_type> particle_group;
    std::shared_ptr<phase_space_type> phase_space;
    std::shared_ptr<potential_type> potential;
    std::shared_ptr<random_type> random;

    brownian_harmonic();
    void test();
};

/** solve the differential equation @f$ \dot r = v = const. @f$ */
template <typename modules_type>
void brownian_harmonic<modules_type>::test()
{
    // construct blocking scheme module (contains the logic)
    observables::dynamics::blocking_scheme blocking_scheme(
        clock, maximum_lag_time, resolution, block_size, 1, 100
    );

    // allocate space for block samples
    auto block_sample = std::make_shared<observables::samples::blocking_scheme<sample_type>>(
        [=]() { return phase_space->template acquire<sample_type>("position"); }
      , blocking_scheme.count()
      , blocking_scheme.block_size()
    );

    // construct correlation functions and connect to blocking scheme logic
    auto correlation_msd = make_shared<observables::dynamics::correlation<msd_type>>(
        make_shared<msd_type>(), block_sample, block_sample
    );
    blocking_scheme.on_correlate(correlation_msd);

    auto correlation_mqd = make_shared<observables::dynamics::correlation<mqd_type>>(
        make_shared<mqd_type>(), block_sample, block_sample
    );
    blocking_scheme.on_correlate(correlation_mqd);

    blocking_scheme.on_sample(block_sample);

    BOOST_TEST_MESSAGE("run Brownian integrator over " << steps << " steps");
    for (size_t i = 0; i < steps; ++i) {
        integrator->integrate();
        clock->advance();
        blocking_scheme.sample();
    }
    blocking_scheme.finalise();

    auto time = blocking_scheme.time()[0];
    auto msd = correlation_msd->result()[0];
    auto mqd = correlation_mqd->result()[0];

    for (size_t i = 0; i < size_t(maximum_lag_time / timestep); ++i) {
        BOOST_TEST_MESSAGE(time[i] << ": " << mean(msd[i]));
        BOOST_CHECK_CLOSE_FRACTION(
            mean(msd[i])
          , 2 * dimension * temperature / stiffness(0) * (1 - expf(-diff_const(0,0) / temperature * stiffness(0) * time[i]))
          , 4.5 * abs(error_of_mean(msd[i]))
        );
    }
}

/**
 * Initialize integrator and dependencies, set basic parameters.
 */
template <typename modules_type>
brownian_harmonic<modules_type>::brownian_harmonic()
{
    timestep = 0.01;
    maximum_lag_time = 7;

    steps = (gpu ? 100 : 100) * maximum_lag_time / timestep;
    resolution = 0.01;
    block_size = 10000;

    temperature = 1;
    npart = gpu ? 10 : 10;

    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = 10.;
    }

    // diffusion constant
    // @param1 displacement
    // @param2 rotational
    diff_const = typename integrator_type::matrix_type(1, 2);
    diff_const <<= 1.0, 1.0;

    stiffness = typename potential_type::scalar_container_type(1);
    stiffness <<= 1.;

    typename potential_type::vector_container_type offset(1);
    offset <<= typename potential_type::vector_type(0);

    particle = std::make_shared<particle_type>(npart, 1);
    particle_group = std::make_shared<particle_group_type>(particle);
    box = std::make_shared<box_type>(edges);
    random = std::make_shared<random_type>();
    potential = std::make_shared<potential_type>(stiffness, offset);
    force = std::make_shared<force_type>(potential, particle, box);
    particle->on_prepend_force([=](){ force->check_cache(); });
    particle->on_force([=](){ force->apply(); });
    integrator = std::make_shared<integrator_type>(particle, random, box, timestep, temperature, diff_const);
    clock = std::make_shared<clock_type>();
    clock->set_timestep(integrator->timestep());
    phase_space = std::make_shared<phase_space_type>(particle, particle_group, box);
}

/**
 * Specify concretely which modules to use: Host modules.
 */
template <int dimension, typename float_type>
struct host_modules_free
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::host::integrators::brownian<dimension, float_type> integrator_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef observables::host::samples::sample<dimension, float_type> sample_type;
    typedef observables::host::phase_space<dimension, float_type> phase_space_type;
    typedef observables::host::dynamics::mean_square_displacement<dimension, float_type> msd_type;
    typedef observables::host::dynamics::mean_quartic_displacement<dimension, float_type> mqd_type;

    static bool const gpu = false;
};

template <int dimension, typename float_type>
struct host_modules_harmonic
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::host::potentials::external::harmonic<dimension, float_type> potential_type;
    typedef mdsim::host::forces::external<dimension, float_type, potential_type> force_type;
    typedef mdsim::host::integrators::brownian<dimension, float_type> integrator_type;
    typedef halmd::random::host::random random_type;
    typedef observables::host::samples::sample<dimension, float_type> sample_type;
    typedef observables::host::phase_space<dimension, float_type> phase_space_type;
    typedef observables::host::dynamics::mean_square_displacement<dimension, float_type> msd_type;
    typedef observables::host::dynamics::mean_quartic_displacement<dimension, float_type> mqd_type;

    static bool const gpu = false;
};

#ifdef HALMD_WITH_GPU

/**
 * Specify concretely which modules to use: GPU modules.
 */
template <int dimension, typename float_type>
struct gpu_modules_free
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef mdsim::gpu::integrators::brownian<dimension, float_type, halmd::random::gpu::rand48> integrator_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef observables::gpu::samples::sample<dimension, float4> sample_type;
    typedef observables::gpu::phase_space<dimension, float_type> phase_space_type;
    typedef observables::gpu::dynamics::mean_square_displacement<dimension, float4> msd_type;
    typedef observables::gpu::dynamics::mean_quartic_displacement<dimension, float4> mqd_type;

    static bool const gpu = true;
};

template <int dimension, typename float_type>
struct gpu_modules_harmonic
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::gpu::potentials::external::harmonic<dimension, float> potential_type;
    typedef mdsim::gpu::forces::external<dimension, float_type, potential_type> force_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef mdsim::gpu::integrators::brownian<dimension, float_type, halmd::random::gpu::rand48> integrator_type;
    typedef observables::gpu::samples::sample<dimension, float4> sample_type;
    typedef observables::gpu::phase_space<dimension, float_type> phase_space_type;
    typedef observables::gpu::dynamics::mean_square_displacement<dimension, float4> msd_type;
    typedef observables::gpu::dynamics::mean_quartic_displacement<dimension, float4> mqd_type;

    static bool const gpu = true;
};

#endif // HALMD_WITH_GPU

#ifndef USE_HOST_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE(brownian_free_host_2d) {
    brownian_free<host_modules_free<2, double>>().test();
}
BOOST_AUTO_TEST_CASE(brownian_free_host_3d) {
    brownian_free<host_modules_free<3, double>>().test();
}
BOOST_AUTO_TEST_CASE(brownian_harmonic_host_2d) {
    brownian_harmonic<host_modules_harmonic<2, double>>().test();
}
BOOST_AUTO_TEST_CASE(brownian_harmonic_host_3d) {
    brownian_harmonic<host_modules_harmonic<3, double>>().test();
}
#else
BOOST_AUTO_TEST_CASE(brownian_free_host_2d) {
    brownian_free<host_modules_free<2, float>>().test();
}
BOOST_AUTO_TEST_CASE(brownian_free_host_3d) {
    brownian_free<host_modules_free<3, float>>().test();
}
BOOST_AUTO_TEST_CASE(brownian_harmonic_host_2d) {
    brownian_harmonic<host_modules_harmonic<2, float>>().test();
}
BOOST_AUTO_TEST_CASE(brownian_harmonic_host_3d) {
    brownian_harmonic<host_modules_harmonic<3, float>>().test();
}
#endif

#ifdef HALMD_WITH_GPU
# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE(brownian_free_gpu_2d, device) {
    brownian_free<gpu_modules_free<2, float>>().test();
}
BOOST_FIXTURE_TEST_CASE(brownian_free_gpu_3d, device) {
    brownian_free<gpu_modules_free<3, float>>().test();
}
BOOST_FIXTURE_TEST_CASE(brownian_harmonic_gpu_2d, device) {
    brownian_harmonic<gpu_modules_harmonic<2, float>>().test();
}
BOOST_FIXTURE_TEST_CASE(brownian_harmonic_gpu_3d, device) {
    brownian_harmonic<gpu_modules_harmonic<3, float>>().test();
}
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE(brownian_free_gpu_2d, device) {
    brownian_free<gpu_modules_free<2, dsfloat>>().test();
}
BOOST_FIXTURE_TEST_CASE(brownian_free_gpu_3d, device) {
    brownian_free<gpu_modules_free<3, dsfloat>>().test();
}
BOOST_FIXTURE_TEST_CASE(brownian_harmonic_gpu_2d, device) {
    brownian_harmonic<gpu_modules_harmonic<2, dsfloat>>().test();
}
BOOST_FIXTURE_TEST_CASE(brownian_harmonic_gpu_3d, device) {
    brownian_harmonic<gpu_modules_harmonic<3, dsfloat>>().test();
}
# endif
#endif // HALMD_WITH_GPU
