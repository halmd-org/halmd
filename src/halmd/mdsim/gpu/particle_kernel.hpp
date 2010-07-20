/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP
#define HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP

#include <cuda_wrapper.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension>
struct particle_wrapper
{
    cuda::symbol<unsigned int> nbox;
    cuda::symbol<unsigned int> ntype;
    cuda::texture<unsigned int> ntypes;
    cuda::function<void (float4*, float4*)> tag;
    static particle_wrapper const kernel;
};

template <int dimension>
particle_wrapper<dimension> const& get_particle_kernel()
{
    return particle_wrapper<dimension>::kernel;
}

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP */
