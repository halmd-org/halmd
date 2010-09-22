/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/integrators/verlet.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost::fusion;

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type>
void verlet<dimension, float_type>::select(po::options const& vm)
{
    if (vm["integrator"].as<std::string>() != "verlet") {
        throw unsuitable_module("mismatching option integrator");
    }
}

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // reference CUDA C++ verlet_wrapper
  , wrapper(&verlet_wrapper<dimension>::wrapper)
  // set parameters
  , timestep_half_(0.5 * timestep_)
{
    /*@{ FIXME remove pre-Lua hack */
    shared_ptr<profiler_type> profiler(modules::fetch<profiler_type>(factory, vm));
    register_runtimes(*profiler);
    /*@}*/

#ifdef USE_VERLET_DSFUN
    //
    // Double-single precision requires two single precision
    // "words" per coordinate. We use the first part of a GPU
    // vector for the higher (most significant) words of all
    // particle positions or velocities, and the second part for
    // the lower (least significant) words.
    //
    // The additional memory is allocated using reserve(), which
    // increases the capacity() without changing the size().
    //
    // Take care to pass capacity() as an argument to cuda::copy
    // or cuda::memset calls if needed, as the lower words will
    // be ignored in the operation.
    //
    LOG("using velocity-Verlet integration in double-single precision");
    particle->g_r.reserve(2 * particle->dim.threads());
    // particle images remain in single precision as they
    // contain integer values (and otherwise would not matter
    // for the long-time stability of the Verlet integrator)
    particle->g_v.reserve(2 * particle->dim.threads());
#else
    LOG_WARNING("using velocity-Verlet integration in single precision");
#endif

    try {
        cuda::copy(timestep_, wrapper->timestep);
        cuda::copy(static_cast<vector_type>(box->length()), wrapper->box_length);
    }
    catch (cuda::error const& e) {
        LOG_ERROR(e.what());
        throw exception("failed to initialize Verlet integrator symbols");
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    try {
        scoped_timer<timer> timer_(at_key<integrate_>(runtime_));
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper->integrate(
            particle->g_r, particle->g_image, particle->g_v, particle->g_f);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to stream first leapfrog step on GPU");
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::finalize()
{
    // TODO: possibly a performance critical issue:
    // the old implementation had this loop included in update_forces(),
    // which saves one additional read of the forces plus the additional kernel execution
    // and scheduling
    try {
        scoped_timer<timer> timer_(at_key<finalize_>(runtime_));
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper->finalize(particle->g_v, particle->g_f);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to stream second leapfrog step on GPU");
    }
}

// explicit instantiation
template class verlet<3, float>;
template class verlet<2, float>;

}}} // namespace mdsim::gpu::integrator

template class module<mdsim::gpu::integrators::verlet<3, float> >;
template class module<mdsim::gpu::integrators::verlet<2, float> >;

} // namespace halmd
