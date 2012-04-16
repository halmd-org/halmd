/*
 * Copyright © 2011  Felix Höfling
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

#define BOOST_TEST_MODULE phase_space
#include <boost/test/unit_test.hpp>

#include <algorithm> // std::max
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#ifdef WITH_CUDA
# include <cuda_wrapper/cuda_wrapper.hpp>
# include <halmd/mdsim/gpu/particle_group.hpp>
# include <halmd/mdsim/gpu/particle_kernel.cuh>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/observables/gpu/samples/phase_space.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace std;

/**
 * test acquisition of phase space samples
 */

/**
 * copy GPU sample to host sample, do nothing if input is host sample
 */
#ifdef WITH_CUDA
template <int dimension, typename float_type>
shared_ptr<observables::host::samples::phase_space<dimension, float_type> const>
copy_sample(shared_ptr<observables::gpu::samples::phase_space<dimension, float_type> const> sample)
{
    using mdsim::gpu::particle_kernel::untagged;

    typedef observables::host::samples::phase_space<dimension, float_type> host_sample_type;
    typedef observables::gpu::samples::phase_space<dimension, float_type> gpu_sample_type;
    typedef typename gpu_sample_type::position_array_type::value_type gpu_vector_type;
    typedef typename host_sample_type::position_array_type::value_type vector_type;

    // allocate memory
    shared_ptr<host_sample_type> result = make_shared<host_sample_type>(sample->position().size());
    cuda::host::vector<gpu_vector_type> h_buf(sample->position().size());

    // copy from GPU to host via page-locked memory

    // positions and types
    cuda::copy(sample->position(), h_buf);
    cuda::thread::synchronize();
    for (size_t i = 0; i < h_buf.size(); ++i) {
        tie(result->position()[i], result->species()[i]) = untagged<vector_type>(h_buf[i]);
    }

    // velocities
    cuda::copy(sample->velocity(), h_buf);
    cuda::thread::synchronize();
    std::copy(h_buf.begin(), h_buf.end(), result->velocity().begin());

    return result;
}
#endif


template <int dimension, typename float_type>
shared_ptr<observables::host::samples::phase_space<dimension, float_type> const>
copy_sample(shared_ptr<observables::host::samples::phase_space<dimension, float_type> const> sample)
{
    return sample;
}

template <typename modules_type>
struct phase_space
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename particle_group_type::particle_type particle_type;
    typedef typename modules_type::input_phase_space_type input_phase_space_type;
    typedef typename input_phase_space_type::sample_type input_sample_type;
    typedef typename modules_type::output_phase_space_type output_phase_space_type;
    typedef typename output_phase_space_type::sample_type output_sample_type;
    static bool const gpu = modules_type::gpu;

    typedef mdsim::clock clock_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    enum { dimension = vector_type::static_size };

    vector<unsigned int> npart;
    typename box_type::vector_type box_length;

    shared_ptr<box_type> box;
    shared_ptr<clock_type> clock;
    shared_ptr<particle_type> particle;
    shared_ptr<input_sample_type> input_sample;

    void test();
    phase_space();
};

template <typename modules_type>
void phase_space<modules_type>::test()
{
    float_type const epsilon = numeric_limits<float_type>::epsilon();

    typename input_sample_type::position_array_type& input_position = input_sample->position();
    typename input_sample_type::velocity_array_type& input_velocity = input_sample->velocity();
    typename input_sample_type::species_array_type& input_species = input_sample->species();

    // prepare input sample
    BOOST_CHECK_EQUAL(input_position.size(), accumulate(npart.begin(), npart.end(), 0u));
    BOOST_CHECK_EQUAL(input_velocity.size(), accumulate(npart.begin(), npart.end(), 0u));
    for (unsigned int i = 0, n = 0; i < npart.size(); ++i) { // iterate over particle species
        for (unsigned int j = 0; j < npart[i]; ++n, ++j) { // iterate over particles
            vector_type& r = input_position[n];
            vector_type& v = input_velocity[n];
            unsigned int& type = input_species[n];
            r[0] = float_type(j) + float_type(1) / (i + 1); //< a large, non-integer value
            r[1] = 0;
            r[dimension - 1] = - static_cast<float_type>(j);
            v[0] = static_cast<float_type>(i);
            v[1] = 0;
            v[dimension - 1] = float_type(1) / (j + 1);
            type = i;
        }
    }

    // copy input sample to particle
    input_phase_space_type(make_shared<particle_group_type>(particle), particle, box, clock).set(input_sample);

    // randomly permute particles in memory
    // TODO

    // acquire sample from particle, construct temporary sampler module
    clock->advance();
    shared_ptr<output_sample_type const> output_sample = output_phase_space_type(make_shared<particle_group_type>(particle), particle, box, clock).acquire();
    BOOST_CHECK(output_sample->step() == 1);

    // compare output and input, copy GPU sample to host before
    shared_ptr<input_sample_type const> result = copy_sample(output_sample);

    typename input_sample_type::position_array_type const& result_position = result->position();
    typename input_sample_type::velocity_array_type const& result_velocity = result->velocity();
    typename input_sample_type::species_array_type const& result_species = result->species();

    BOOST_CHECK_EQUAL(result_position.size(), accumulate(npart.begin(), npart.end(), 0u));
    for (unsigned int i = 0, n = 0; i < npart.size(); ++i) { // iterate over particle species
        for (unsigned int j = 0; j < npart[i]; ++n, ++j) { // iterate over particles
            // compare positions with a tolerance due to mapping to and from the periodic box
            for (unsigned int k = 0; k < dimension; ++k) {
                BOOST_CHECK_CLOSE_FRACTION(result_position[n][k], input_position[n][k], 10 * epsilon);
            }
        }
    }
    // compare velocities directly as they should not have been modified
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result_velocity.begin(), result_velocity.end()
      , input_velocity.begin(), input_velocity.end()
    );
    // compare particle species
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result_species.begin(), result_species.end()
      , input_species.begin(), input_species.end()
    );
}

template <typename modules_type>
phase_space<modules_type>::phase_space()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    npart = list_of(1024)(512)(30)(1); //< choose a value smaller than warp size and some limiting values
    vector<double> mass = list_of(1)(1)(1)(1);

    // choose a box length with is not an exactly representable as a
    // floating-point number and which is small enough to have some overflow
    // from the periodic box. In addition, some of the coordinates should sit
    // precisely at the edge.
    box_length = fixed_vector<double, dimension>(40./3);

    // create modules
    particle = make_shared<particle_type>(npart, mass);
    box = make_shared<box_type>(particle->nbox, box_length);
    input_sample = make_shared<input_sample_type>(particle->nbox);
    clock = make_shared<clock_type>(0); // bogus time-step

    // set particle tags and types
    particle->set();
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle_group_all<dimension, float_type> particle_group_type;
    typedef observables::host::phase_space<dimension, float_type> input_phase_space_type;
    typedef input_phase_space_type output_phase_space_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( phase_space_host_2d ) {
    phase_space<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( phase_space_host_3d ) {
    phase_space<host_modules<3, double> >().test();
}

#ifdef WITH_CUDA
template <int dimension, typename float_type>
struct gpu_host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle_group_all<dimension, float_type> particle_group_type;
    typedef observables::gpu::phase_space<observables::host::samples::phase_space<dimension, float_type> > input_phase_space_type;
    typedef input_phase_space_type output_phase_space_type;
    static bool const gpu = true;
};

template <int dimension, typename float_type>
struct gpu_gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle_group_all<dimension, float_type> particle_group_type;
    typedef observables::gpu::phase_space<observables::host::samples::phase_space<dimension, float_type> > input_phase_space_type;
    typedef observables::gpu::phase_space<observables::gpu::samples::phase_space<dimension, float_type> > output_phase_space_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( phase_space_gpu_host_2d, device ) {
    phase_space<gpu_host_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( phase_space_gpu_host_3d, device ) {
    phase_space<gpu_host_modules<3, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( phase_space_gpu_gpu_2d, device ) {
    phase_space<gpu_gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( phase_space_gpu_gpu_3d, device ) {
    phase_space<gpu_gpu_modules<3, float> >().test();
}
#endif // WITH_CUDA
