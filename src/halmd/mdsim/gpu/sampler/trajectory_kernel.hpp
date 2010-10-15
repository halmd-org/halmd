/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_SAMPLE_TRAJECTORY_WRAPPER_CUH
#define HALMD_MDSIM_GPU_SAMPLE_TRAJECTORY_WRAPPER_CUH

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace sampler
{

template <int dimension>
struct trajectory_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    /** positions, types */
    cuda::texture<float4> r;
    /** minimum image vectors */
    cuda::texture<coalesced_vector_type> image;
    /** velocities, tags */
    cuda::texture<float4> v;
    /** cubic box edgle length */
    cuda::symbol<vector_type> box_length;
    /** sample trajectory for all particle of a single species */
    cuda::function<void (unsigned int const*, coalesced_vector_type*, coalesced_vector_type*)> sample;

    static trajectory_wrapper const kernel;
};

template <int dimension>
trajectory_wrapper<dimension> const& get_trajectory_kernel()
{
    return trajectory_wrapper<dimension>::kernel;
}

}}} // namespace mdsim::gpu::sampler

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_SAMPLE_TRAJECTORY_WRAPPER_CUH */
