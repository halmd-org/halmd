/*
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
#include <halmd/mdsim/host/integrators/brownian.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
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
# include <halmd/mdsim/gpu/integrators/brownian.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/orientations/uniform.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/observables/gpu/dynamics/orientational_autocorrelation.hpp>
# include <halmd/observables/gpu/dynamics/mean_quartic_displacement.hpp>
# include <halmd/observables/gpu/dynamics/mean_square_displacement.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace halmd;
using namespace std;

template <typename modules_type>
struct test_brownian
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::integrator_type integrator_type;
    typedef typename integrator_type::matrix_type matrix_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::orientation_type orientation_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::velocity_type velocity_type;
    typedef mdsim::clock clock_type;
    typedef typename modules_type::sample_type sample_type;
    typedef typename modules_type::phase_space_type phase_space_type;

    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;

    typedef typename modules_type::numeric_limits numeric_limits;
    static bool const gpu = modules_type::gpu;

    size_t steps;
    double density;
    double temp;
    double timestep;
    matrix_type D = matrix_type(4, 1);
    double maximum_lag_time;
    double resolution;
    unsigned int block_size;
    unsigned int npart;
    fixed_vector<double, dimension> box_ratios;
    fixed_vector<double, dimension> slab;

    std::shared_ptr<box_type> box;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<integrator_type> integrator;
    std::shared_ptr<random_type> random;
    std::shared_ptr<position_type> position;
    std::shared_ptr<orientation_type> orientation;
    std::shared_ptr<velocity_type> velocity;
    std::shared_ptr<clock_type> clock;
    std::shared_ptr<phase_space_type> phase_space;

    typedef observables::gpu::dynamics::mean_square_displacement<dimension, float_type> msd_type;
    typedef observables::gpu::dynamics::mean_quartic_displacement<dimension, float_type> mqd_type;
    typedef observables::gpu::dynamics::orientational_autocorrelation<dimension, float_type> ocf_type;

    test_brownian();
    void free_brownian_motion();
};

/** solve the differential equation @f$ \dot r = v = const @f$ */
template <typename modules_type>
void test_brownian<modules_type>::free_brownian_motion()
{
    // copy initial positions and velocities from particle to host sample
    std::shared_ptr<sample_type const> initial_sample = phase_space->acquire();

    // construct blocking scheme module (contains the logic)
    observables::dynamics::blocking_scheme blocking_scheme(
        clock, maximum_lag_time, resolution, block_size, 1
    );

    // allocate space for block samples
    auto block_sample = std::make_shared<observables::samples::blocking_scheme<sample_type>>(
        [=]() { return phase_space->acquire(); }, 100, block_size
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

    auto correlation_ocf = make_shared<observables::dynamics::correlation<ocf_type>>(
        make_shared<ocf_type>(), block_sample, block_sample
    );
    blocking_scheme.on_correlate(correlation_ocf);

    blocking_scheme.on_sample(block_sample);
    // perform integration
    BOOST_TEST_MESSAGE("run Brownian integration for free motion over " << steps << " steps");
    if (modules_type::gpu) {
        BOOST_TEST_MESSAGE("number of blocks:  " << random->blocks() << "number of threads:  " << random->threads() );
    }
    for (size_t i = 0; i < steps; ++i) {
        integrator->integrate();
        clock->advance();
        blocking_scheme.sample();
    }
    blocking_scheme.finalise(); 
    BOOST_TEST_MESSAGE("test correlation function of the motion");

    auto time = blocking_scheme.time()[0];
    auto msd = correlation_msd->result()[0];
    auto mqd = correlation_mqd->result()[0];
    auto ocf = correlation_ocf->result()[0];
    
    BOOST_TEST_MESSAGE("check if msd is close to 6 D t");
    double max_deviation = 0;
    for (size_t i = 0; i < time.size(); ++i) {
      //BOOST_TEST_MESSAGE(abs(mean(msd[i])- 6 * time[i]) - 6 * error_of_mean(msd[i])  );
      //BOOST_CHECK_CLOSE(mean(msd[i]), 6 * time[i], 6 * error_of_mean(msd[i]));
      //max_deviation = std::max(abs(mean(msd[i]) - 6 * time[i]), max_deviation);
      BOOST_TEST_MESSAGE( time[i] << " " << mean(msd[i]) <<  "  " << mean(msd[i]) / time[i] << "  " << 6 * time[i] << "  " << error_of_mean(msd[i]) << "  " << mean(ocf[i]) << "  " << (3 * mean(mqd[i])   / (5 * mean(msd[i]) * mean(msd[i]))-1) ); // << "  " << (3 * error_of_mean(mqd[i])   / (5 * pow(mean(msd[i]), 2))) + (6 * mean(mqd[i]) * error_of_mean(mqd[i])  / (5 * pow(mean(msd[i]), 3))));
    }
    BOOST_TEST_MESSAGE("Maximum deviation  " << max_deviation);
}

/**
 * Initialize integrator and dependencies, set basic parameters.
 */
template <typename modules_type>
test_brownian<modules_type>::test_brownian()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");
    typedef fixed_vector<double, dimension> vector_type;

    // set test parameters
    steps = 10000; // run for as many steps as possible, wrap around the box for about 10 times
    maximum_lag_time = 10;
    resolution = 0.01;
    block_size = 1000;
    // set module parameters
    density = 1e-6; // a low density implies large values of the position vectors
    temp = 1; // the temperature defines the average velocities
    timestep = 0.001; // small, but typical timestep
    npart = gpu ? 512 : 108; // optimize filling of fcc lattice, use only few particles on the host
    box_ratios = (dimension == 3) ? vector_type{1., 2., 1.01} : vector_type{1., 2.};
    //unsigned int seed = 1; 
    double det = accumulate(box_ratios.begin(), box_ratios.end(), 1., multiplies<double>());
    double volume = npart / density;
    double edge_length = pow(volume / det, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = edge_length * box_ratios[i];
    }
    slab = 1;
    D <<= 1.0, 2.0, 1.0, 30.0;
    // create modules
    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    random = std::make_shared<random_type>(365324873, 2, 512);
    integrator = std::make_shared<integrator_type>(particle, random, box, timestep, temp, D);
    position = std::make_shared<position_type>(particle, box, slab);
    orientation = std::make_shared<orientation_type>(particle, random);
    velocity = std::make_shared<velocity_type>(particle, random, temp);
    clock = std::make_shared<clock_type>();
    clock->set_timestep(integrator->timestep());
    std::shared_ptr<particle_group_type> particle_group = std::make_shared<particle_group_type>(particle);
    phase_space = std::make_shared<phase_space_type>(particle, particle_group, box, clock);

    // set positions and velocities
    BOOST_TEST_MESSAGE("position particles on lattice");
    position->set();
    BOOST_TEST_MESSAGE("position orientation uniformly");
    orientation->set();
    BOOST_TEST_MESSAGE("set particle velocities");
    velocity->set();
    BOOST_TEST_MESSAGE("using timestep  " << integrator->timestep());
}

/**
 * Specify concretely which modules to use: Host modules.
 */
template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::host::integrators::brownian<dimension, float_type> integrator_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
    typedef observables::host::phase_space<dimension, float_type> phase_space_type;

    typedef typename std::numeric_limits<float_type> numeric_limits;

    static bool const gpu = false;

    static void set_velocity(std::shared_ptr<particle_type> particle);
};

/** host specific helper function: set particle velocity to v = -r */
template <int dimension, typename float_type>
void host_modules<dimension, float_type>::set_velocity(std::shared_ptr<particle_type> particle)
{
    // copy -r[i] to v[i]
    auto const& position = read_cache(particle->position());
    auto velocity = make_cache_mutable(particle->velocity());
    std::transform(
        position.begin(), position.end()
      , velocity->begin()
      , negate<typename particle_type::vector_type>()
    );
}
#ifdef HALMD_WITH_GPU
/**
 * Specify concretely which modules to use: Gpu modules.
 */
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::mrg32k3a> random_type;
    typedef mdsim::gpu::integrators::brownian<dimension, float_type, halmd::random::gpu::mrg32k3a> integrator_type;
    typedef mdsim::gpu::orientations::uniform<dimension, float_type, halmd::random::gpu::mrg32k3a> orientation_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::mrg32k3a> velocity_type;
    //typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
    typedef observables::gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef observables::gpu::phase_space<sample_type> phase_space_type;

    static bool const gpu = true;

#ifndef USE_VERLET_DSFUN
    typedef typename std::numeric_limits<float_type> numeric_limits;
#else
    // FIXME define numeric_limits for dsfloat
    // see, e.g., http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
    struct numeric_limits {
        static float_type epsilon() { return std::pow(float_type(2), -44); }
        static float_type min() { return std::numeric_limits<float>::min(); }
    };
#endif

    static void set_velocity(std::shared_ptr<particle_type> particle);
};

/** GPU specific helper function: set particle velocity to v = -r */
template <int dimension, typename float_type>
void gpu_modules<dimension, float_type>::set_velocity(std::shared_ptr<particle_type> particle)
{
    using namespace algorithm::gpu;
    typedef typename particle_type::vector_type vector_type;
    typedef apply_wrapper<negate_, vector_type, float4, vector_type, float4> apply_negate_wrapper;

    auto const& position = read_cache(particle->position());
    auto velocity = make_cache_mutable(particle->velocity());

    // copy -g_r[i] to g_v[i]
    //
    // The kernel wrapper is declared in algorithm/gpu/apply_kernel.hpp
    // and the CUDA kernel is instantiated in apply_negate.cu.
    //
    // If the arrays are stored as two subsequent float4 arrays for
    // double-single representation, we the negation is applied to both floats
    // independently (in correspondence to the definition of operator- for
    // dsfloat).
    //
    // Caveat: overwrites particle tags in g_v (which are not used anyway)
    try {
        cuda::configure(particle->dim.grid, particle->dim.block);
        apply_negate_wrapper::kernel.apply(&*position.begin(), &*velocity->begin(), position.capacity());
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("copying negated positions to velocities on GPU failed");
        throw;
    }
}

BOOST_FIXTURE_TEST_CASE( brownian_gpu_3d_overdamped, device ) {
    test_brownian<gpu_modules<3, float> >().free_brownian_motion();
}
#endif // HALMD_WITH_GPU
/*
BOOST_AUTO_TEST_CASE( brownian_host_3d_overdamped ) {
    test_brownian<host_modules<3, double> >().free_brownian_motion();
}
*/
