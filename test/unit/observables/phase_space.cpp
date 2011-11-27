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
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/phase_space.hpp>
#include <halmd/mdsim/host/velocities/phase_space.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#ifdef WITH_CUDA
# include <cuda_wrapper/cuda_wrapper.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_kernel.cuh>
# include <halmd/mdsim/gpu/positions/phase_space.hpp>
# include <halmd/mdsim/gpu/velocities/phase_space.hpp>
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
shared_ptr<observables::host::samples::phase_space<dimension, float_type> >
copy_sample(shared_ptr<observables::gpu::samples::phase_space<dimension, float_type> > sample)
{
    using mdsim::gpu::particle_kernel::untagged;

    typedef observables::host::samples::phase_space<dimension, float_type> host_sample_type;
    typedef typename observables::gpu::samples::phase_space<dimension, float_type>::gpu_vector_type gpu_vector_type;
    typedef typename observables::gpu::samples::phase_space<dimension, float_type>::vector_type vector_type;

    // allocate memory
    shared_ptr<host_sample_type> result = make_shared<host_sample_type>(sample->r->size());
    cuda::host::vector<gpu_vector_type> h_buf(sample->r->size());

    // copy from GPU to host via page-locked memory

    // positions and types
    cuda::copy(*sample->r, h_buf);
    cuda::thread::synchronize();
    for (size_t i = 0; i < h_buf.size(); ++i) {
        tie((*result->r)[i], (*result->type)[i]) = untagged<vector_type>(h_buf[i]);
    }

    // velocities
    cuda::copy(*sample->v, h_buf);
    cuda::thread::synchronize();
    std::copy(h_buf.begin(), h_buf.end(), result->v->begin());

    return result;
}
#endif


template <int dimension, typename float_type>
shared_ptr<observables::host::samples::phase_space<dimension, float_type> >
copy_sample(shared_ptr<observables::host::samples::phase_space<dimension, float_type> > sample)
{
    return sample;
}

template <typename modules_type>
struct phase_space
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::velocity_type velocity_type;
    typedef typename modules_type::phase_space_type phase_space_type;
    typedef typename modules_type::input_sample_type input_sample_type;
    typedef typename modules_type::output_sample_type output_sample_type;
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
    shared_ptr<position_type> position;
    shared_ptr<velocity_type> velocity;
    shared_ptr<input_sample_type> input_sample;
    shared_ptr<output_sample_type> output_sample;

    void test();
    phase_space();
};

template <typename modules_type>
void phase_space<modules_type>::test()
{
    float_type const epsilon = numeric_limits<float_type>::epsilon();

    // prepare input sample
    BOOST_CHECK_EQUAL(input_sample->r->size(), accumulate(npart.begin(), npart.end(), 0));
    BOOST_CHECK_EQUAL(input_sample->v->size(), accumulate(npart.begin(), npart.end(), 0));
    for (unsigned int i = 0, n = 0; i < npart.size(); ++i) { // iterate over particle species
        for (unsigned int j = 0; j < npart[i]; ++n, ++j) { // iterate over particles
            vector_type& r = (*input_sample->r)[n];
            vector_type& v = (*input_sample->v)[n];
            unsigned int& type = (*input_sample->type)[n];
            r[0] = float_type(j) + float_type(1) / (i + 1); //< a large, non-integer value
            r[1] = 0;
            r[dimension - 1] = - static_cast<float_type>(j);
            v[0] = static_cast<float_type>(i);
            v[1] = 0;
            v[dimension - 1] = float_type(1) / (j + 1);
            type = i;
        }
    }
    input_sample->step = 0;

    // copy input sample to particle
    position->set();
    velocity->set();

    // randomly permute particles in memory
    // TODO

    // acquire sample from particle, construct temporary sampler module
    clock->advance();
    phase_space_type(output_sample, particle, box, clock).acquire();
    BOOST_CHECK(output_sample->step == 1);

    // compare output and input, copy GPU sample to host before
    shared_ptr<observables::host::samples::phase_space<dimension, float_type> > result
        = copy_sample(output_sample);
    BOOST_CHECK_EQUAL(result->r->size(), accumulate(npart.begin(), npart.end(), 0));
    for (unsigned int i = 0, n = 0; i < npart.size(); ++i) { // iterate over particle species
        for (unsigned int j = 0; j < npart[i]; ++n, ++j) { // iterate over particles
            // compare positions with a tolerance due to mapping to and from the periodic box
            for (unsigned int k = 0; k < dimension; ++k) {
                BOOST_CHECK_CLOSE_FRACTION((*result->r)[n][k], (*input_sample->r)[n][k], 10 * epsilon);
            }
        }
    }
    // compare velocities directly as they should not have been modified
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result->v->begin(), result->v->end()
      , input_sample->v->begin(), input_sample->v->end()
    );
    // compare particle species
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result->type->begin(), result->type->end()
      , input_sample->type->begin(), input_sample->type->end()
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
    output_sample = make_shared<output_sample_type>(particle->nbox);
    position = make_shared<position_type>(particle, box, input_sample);
    velocity = make_shared<velocity_type>(particle, input_sample);
    clock = make_shared<clock_type>(0); // bogus time-step

    // set particle tags and types
    particle->set();
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::positions::phase_space<dimension, float_type> position_type;
    typedef mdsim::host::velocities::phase_space<dimension, float_type> velocity_type;
    typedef observables::host::phase_space<dimension, float_type> phase_space_type;
    typedef observables::host::samples::phase_space<dimension, float_type> input_sample_type;
    typedef input_sample_type output_sample_type;
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
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::positions::phase_space<dimension, float_type> position_type;
    typedef mdsim::gpu::velocities::phase_space<dimension, float_type> velocity_type;
    typedef observables::host::samples::phase_space<dimension, float_type> input_sample_type;
    typedef input_sample_type output_sample_type;
    typedef observables::gpu::phase_space<output_sample_type> phase_space_type;
    static bool const gpu = true;
};

template <int dimension, typename float_type>
struct gpu_gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::positions::phase_space<dimension, float_type> position_type;
    typedef mdsim::gpu::velocities::phase_space<dimension, float_type> velocity_type;
    typedef observables::host::samples::phase_space<dimension, float_type> input_sample_type;
    typedef observables::gpu::samples::phase_space<dimension, float_type> output_sample_type;
    typedef observables::gpu::phase_space<output_sample_type> phase_space_type;
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
