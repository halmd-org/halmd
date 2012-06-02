/*
 * Copyright © 2011  Michael Kopp and Felix Höfling
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE euler
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/assign.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <functional>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/integrators/euler.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/particle_groups/from_range.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
# include <halmd/algorithm/gpu/apply_kernel.hpp>
# include <halmd/mdsim/gpu/integrators/euler.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace std;

/** test Euler integrator: ordinary differential equations of 1st order
 *
 * Two differential equations are solved using the Euler integrator,
 * @f\[
 *   \dot r = v = \mathit{const} \quad \text{and} \quad  \dot r = -r \, .
 * @f\]
 * The results are then compared to the algebraic solutions of the
 * numerical scheme to properly account for discretisation errors,
 * @f\[
 *   r(n) = r_0 + v \, dt \, n
 *   \quad \text{and} \quad
 *   r(n) = (1 - dt)^n \, r_0 ,
 * @f\]
 * with \f$n\f$ the number of steps taken.
 *
 */

template <typename modules_type>
struct test_euler
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::integrator_type integrator_type;
    typedef typename modules_type::particle_type particle_type;
    typedef mdsim::particle_groups::from_range<particle_type> particle_group_type;
    typedef typename modules_type::position_type position_type;
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
    unsigned int npart;
    fixed_vector<double, dimension> box_ratios;
    fixed_vector<double, dimension> slab;

    std::shared_ptr<box_type> box;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<integrator_type> integrator;
    std::shared_ptr<random_type> random;
    std::shared_ptr<position_type> position;
    std::shared_ptr<velocity_type> velocity;
    std::shared_ptr<clock_type> clock;
    std::shared_ptr<phase_space_type> phase_space;

    test_euler();
    void linear_motion();
    void overdamped_motion();
};

/** solve the differential equation @f$ \dot r = v = const @f$ */
template <typename modules_type>
void test_euler<modules_type>::linear_motion()
{
    // copy initial positions and velocities from particle to host sample
    std::shared_ptr<sample_type const> initial_sample = phase_space->acquire();

    // perform integration
    BOOST_TEST_MESSAGE("running Euler integration for linear motion over " << steps << " steps");
    clock->set_timestep(integrator->timestep());
    for (size_t i = 0; i < steps; ++i) {
        integrator->integrate();
        clock->advance();
    }

    // acquire sample with final positions and velocities
    std::shared_ptr<sample_type const> sample = phase_space->acquire();

    typename sample_type::position_array_type const& initial_position = initial_sample->position();
    typename sample_type::velocity_array_type const& initial_velocity = initial_sample->velocity();
    typename sample_type::position_array_type const& position = sample->position();

    // particlewise comparison with analytic solution
    // the absolute error should be relative to the maximum value, i.e., the box length
    float_type tolerance = 4 * steps * numeric_limits::epsilon() * norm_inf(box->length());
    float_type duration = steps * integrator->timestep();
    float_type max_deviation = 0;
    for (size_t i = 0; i < npart; ++i) {
        vector_type const& r0 = initial_position[i];
        vector_type const& v0 = initial_velocity[i];
        vector_type const& r_final = position[i];

        vector_type r_analytic = r0 + duration * v0;

        // check that maximum deviation of a vector component is less than the tolerance
        BOOST_CHECK_SMALL(norm_inf(r_final - r_analytic), tolerance);
        max_deviation = std::max(norm_inf(r_final - r_analytic), max_deviation);
    }
    BOOST_TEST_MESSAGE("Maximum deviation: " << max_deviation << ", tolerance: " << tolerance);
    BOOST_CHECK_LT(max_deviation, tolerance);
}

/** solve the differential equation @f$ \dot r = - r @f$ */
template <typename modules_type>
void test_euler<modules_type>::overdamped_motion()
{
    // copy initial positions and velocities from particle to host sample
    std::shared_ptr<sample_type const> initial_sample = phase_space->acquire();

    // reduce number of steps as the test runs much slower
    // and the outcome can't be well represented by float
    steps /= gpu ? 100 : 10;

    // perform integration
    BOOST_TEST_MESSAGE("running Euler integration for overdamped motion over " << steps << " steps");
    clock->set_timestep(integrator->timestep());
    for (size_t i = 0; i < steps; ++i) {
        modules_type::set_velocity(particle); // set particle velocity: v = -r
        integrator->integrate();
        clock->advance();
    }

    // acquire sample with final positions and velocities
    std::shared_ptr<sample_type const> sample = phase_space->acquire();

    // particlewise comparison with analytic solution
    // r_n = r_0 * (1 - Δt)^n → r_0 * exp(-n Δt)
    float_type factor = pow(1 - integrator->timestep(), static_cast<double>(steps));
    float_type max_deviation = 0;
    for (size_t i = 0; i < npart; ++i) {
        vector_type const& r0 = initial_sample->position()[i];
        vector_type const& r_final = sample->position()[i];

        vector_type r_analytic = r0 * factor;

        // Check that maximum deviation of a vector component is less than the tolerance
        // the tolerance is computed by summing up all errors:
        // @f$$ E_{total} = \epsilon \, \sum_n x_n = \epsilon (x_0 - x_n) \frac{1-\Delta t}{\Delta t} @f$$
        // and @f$\epsilon @f$ is the relative error for one addition.
        float_type tolerance_step = numeric_limits::epsilon() * norm_inf(r0 - r_final) * (1 - timestep) / timestep;
        tolerance_step = std::max(tolerance_step, numeric_limits::min()); // avoid "0 < 0"
        BOOST_CHECK_SMALL(norm_inf(r_final - r_analytic), tolerance_step);
        max_deviation = std::max(norm_inf(r_final - r_analytic), max_deviation);
    }
    BOOST_TEST_MESSAGE("Maximum deviation: " << max_deviation);
}

/**
 * Initialize integrator and dependencies, set basic parameters.
 */
template <typename modules_type>
test_euler<modules_type>::test_euler()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set test parameters
    steps = 1000000; // run for as many steps as possible, wrap around the box for about 10 times

    // set module parameters
    density = 1e-6; // a low density implies large values of the position vectors
    temp = 1; // the temperature defines the average velocities
    timestep = 0.001; // small, but typical timestep
    npart = gpu ? 4000 : 108; // optimize filling of fcc lattice, use only few particles on the host
    box_ratios = (dimension == 3) ? list_of<double>(1)(2)(1.01) : list_of<double>(1)(2);
    double det = accumulate(box_ratios.begin(), box_ratios.end(), 1., multiplies<double>());
    double volume = npart / density;
    double edge_length = pow(volume / det, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = edge_length * box_ratios[i];
    }
    slab = 1;

    // create modules
    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    integrator = std::make_shared<integrator_type>(particle, box, timestep);
    random = std::make_shared<random_type>();
    position = std::make_shared<position_type>(particle, box, slab);
    velocity = std::make_shared<velocity_type>(particle, random, temp);
    clock = std::make_shared<clock_type>();
    std::shared_ptr<particle_group_type> particle_group = std::make_shared<particle_group_type>(particle, 0, particle->nparticle());
    phase_space = std::make_shared<phase_space_type>(particle, particle_group, box, clock);

    // set positions and velocities
    BOOST_TEST_MESSAGE("position particles on lattice");
    position->set();
    BOOST_TEST_MESSAGE("set particle velocities");
    velocity->set();
}

/**
 * Specify concretely which modules to use: Host modules.
 */
template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::integrators::euler<dimension, float_type> integrator_type;
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
    transform(
        particle->position().begin(), particle->position().end()
      , particle->velocity().begin()
      , negate<typename particle_type::vector_type>()
    );
}

BOOST_AUTO_TEST_CASE( euler_host_2d_linear ) {
    test_euler<host_modules<2, double> >().linear_motion();
}
BOOST_AUTO_TEST_CASE( euler_host_3d_linear ) {
    test_euler<host_modules<3, double> >().linear_motion();
}

BOOST_AUTO_TEST_CASE( euler_host_2d_overdamped ) {
    test_euler<host_modules<2, double> >().overdamped_motion();
}
BOOST_AUTO_TEST_CASE( euler_host_3d_overdamped ) {
    test_euler<host_modules<3, double> >().overdamped_motion();
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
    typedef mdsim::gpu::integrators::euler<dimension, float_type> integrator_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
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
        apply_negate_wrapper::kernel.apply(particle->position(), particle->velocity(), particle->position().capacity());
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("copying negated positions to velocities on GPU failed");
        throw;
    }
}

BOOST_FIXTURE_TEST_CASE( euler_gpu_2d_linear, device ) {
    test_euler<gpu_modules<2, float> >().linear_motion();
}
BOOST_FIXTURE_TEST_CASE( euler_gpu_3d_linear, device ) {
    test_euler<gpu_modules<3, float> >().linear_motion();
}

BOOST_FIXTURE_TEST_CASE( euler_gpu_2d_overdamped, device ) {
    test_euler<gpu_modules<2, float> >().overdamped_motion();
}
BOOST_FIXTURE_TEST_CASE( euler_gpu_3d_overdamped, device ) {
    test_euler<gpu_modules<3, float> >().overdamped_motion();
}
#endif // HALMD_WITH_GPU
