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

#ifndef HALMD_MDSIM_GPU_INTEGRATOR_BROWNIAN_KERNEL_HPP
#define HALMD_MDSIM_GPU_INTEGRATOR_BROWNIAN_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename rng_type>
struct brownian_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_pseudo_vector_type coalesced_pseudo_vector_type;
    typedef typename type_traits<4, float_type>::gpu::ptr_type ptr_type;
    typedef typename type_traits<4, float_type>::gpu::const_ptr_type const_ptr_type;

    static cuda::texture<float2> param;

    typedef cuda::function <void (
        ptr_type
      , coalesced_vector_type*
      , const_ptr_type
      , coalesced_vector_type const*
      , coalesced_pseudo_vector_type const*
      , float
      , float
      , rng_type
      , unsigned int
      , unsigned int
      , vector_type
    )> integrate_kernel_type;

    integrate_kernel_type integrate_position;
    integrate_kernel_type integrate_orientation;
    integrate_kernel_type integrate_both;

    static brownian_wrapper const kernel;
};

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATOR_BROWNIAN_KERNEL_HPP */
