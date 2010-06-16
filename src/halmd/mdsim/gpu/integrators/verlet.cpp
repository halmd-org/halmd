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
#include <halmd/mdsim/gpu/integrators/verlet_kernel.cuh>
#include <halmd/utility/module.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::options(po::options_description& desc)
{
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::depends()
{
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::select(po::options const& vm)
{
    if (vm["integrator"].as<std::string>() != "verlet") {
        throw unsuitable_module("mismatching option integrator");
    }
}

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(po::options const& vm)
  : _Base(vm)
  // set parameters
  , timestep_half_(0.5 * timestep_)
{
    LOG("using velocity-Verlet integration");

    unsigned blocks = (particle->nbox + device->threads() - 1) / device->threads();
    dim_ = cuda::config(blocks, device->threads());
    LOG_DEBUG("number of CUDA execution blocks: " << dim_.blocks_per_grid());
    LOG_DEBUG("number of CUDA execution threads per block: " << dim_.threads_per_block());
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    try {
        cuda::configure(dim_.grid, dim_.block);
        verlet_wrapper<dimension>::integrate(
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
        cuda::configure(dim_.grid, dim_.block);
        verlet_wrapper<dimension>::finalize(particle->g_v, particle->g_f);
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
