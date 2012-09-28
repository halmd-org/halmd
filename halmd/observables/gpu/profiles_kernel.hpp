/*
 * Copyright © 2008-2010  Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_PROFILES_KERNEL_HPP
#define HALMD_OBSERVABLES_GPU_PROFILES_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd
{
namespace observables { namespace gpu
{

template <int dimension>
struct profiles_wrapper
{
    typedef typename mdsim::type_traits<dimension, float> type_traits;
    typedef typename type_traits::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename type_traits::gpu::vector_type gpu_vector_type;
    typedef typename type_traits::vector_type vector_type;
    typedef typename type_traits::gpu::stress_tensor_first_type gpu_stress_tensor_first_type;
    typedef typename mdsim::type_traits<dimension, unsigned>::vector_type index_type;
    typedef typename mdsim::type_traits<dimension, unsigned>::gpu::vector_type gpu_index_type;

    /** number of particles */
    cuda::symbol<unsigned> nbox;

    cuda::function<void (float4 const*, uint*, uint*, index_type, vector_type, vector_type)> compute_bins;
    cuda::function<void (uint const*, int2*)> find_boundaries;
    cuda::function<void (
        int2 const*
      , uint const*
      , gpu_stress_tensor_first_type const*
      , float4 const*
      , coalesced_vector_type*
      , uint
    )> collect_stress_tensor;

    static profiles_wrapper const kernel;
};

template <int dimension>
profiles_wrapper<dimension> const& get_profiles_kernel()
{
    return profiles_wrapper<dimension>::kernel;
}

}} // namespace observables::gpu

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PROFILES_KERNEL_HPP */
