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
#include <boost/algorithm/string/predicate.hpp>
#include <exception>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_kernel.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Resolve module dependencies
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::depends()
{
    // The gpu::device module selects the GPU, which has to occur
    // before the first CUDA kernel call or memory allocation. As
    // the particle module is instantiated foremost and allocates
    // GPU memory, we add a dependency on gpu::device.
    modules::depends<_Self, device_type>::required();
}

template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::select(po::variables_map const& vm)
{
    if (!starts_with(vm["backend"].as<string>(), "gpu")) {
        throw unsuitable_module("mismatching option backend");
    }
}

template <unsigned int dimension, typename float_type>
particle<dimension, float_type>::particle(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // dependency injection
  , device(modules::fetch<device_type>(factory, vm))
  // default CUDA kernel execution dimensions
  , dim(cuda::config((nbox + device->threads() - 1) / device->threads(), device->threads()))
  // allocate global device memory
  , g_r(nbox)
  , g_image(nbox)
  , g_v(nbox)
  , g_f(nbox)
  // allocate page-locked host memory
  , h_r(nbox)
  , h_image(nbox)
  , h_v(nbox)
{
    LOG_DEBUG("number of CUDA execution blocks: " << dim.blocks_per_grid());
    LOG_DEBUG("number of CUDA execution threads per block: " << dim.threads_per_block());

    //
    // As the number of threads may exceed the nmber of particles
    // to account for an integer number of threads per block,
    // we need to allocate excess memory for the GPU vectors.
    //
    // The additional memory is allocated using reserve(), which
    // increases the capacity() without changing the size(). The
    // coordinates of these "virtual" particles will be ignored
    // in cuda::copy or cuda::memset calls.
    //
    try {
        g_r.reserve(dim.threads());
        g_image.reserve(dim.threads());
        g_v.reserve(dim.threads());
        g_f.reserve(dim.threads());
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to allocate particles in global device memory");
    }

    try {
        cuda::copy(nbox, get_particle_kernel<dimension>().nbox);
        cuda::copy(ntype, get_particle_kernel<dimension>().ntype);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to copy particle parameters to device symbols");
    }
}

/**
 * set particle tags and types
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::set()
{
    try {
        cuda::configure(dim.grid, dim.block);
        cuda::vector<unsigned int> g_ntypes(ntypes.size());
        cuda::copy(ntypes, g_ntypes);
        get_particle_kernel<dimension>().ntypes.bind(g_ntypes);
        get_particle_kernel<dimension>().tag(g_r, g_v);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to set particle tags and types");
    }
}

// explicit instantiation
template class particle<3, float>;
template class particle<2, float>;

}} // namespace mdsim::gpu

template class module<mdsim::gpu::particle<3, float> >;
template class module<mdsim::gpu::particle<2, float> >;

} // namespace halmd
