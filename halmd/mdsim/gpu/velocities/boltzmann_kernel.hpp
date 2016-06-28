/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_VELOCITIES_BOLTZMANN_KERNEL_HPP
#define HALMD_MDSIM_GPU_VELOCITIES_BOLTZMANN_KERNEL_HPP

#include <boost/mpl/if.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

//
// Maxwell-Boltzmann distribution at accurate temperature
//

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <
    int dimension
  , typename float_type
  , typename rng_type
>
struct boltzmann_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef cuda::function<void (float4*, unsigned int, unsigned int, float, coalesced_vector_type*, dsfloat*, dsfloat*, rng_type)> gaussian_impl_type;
    gaussian_impl_type gaussian_impl_32;
    gaussian_impl_type gaussian_impl_64;
    gaussian_impl_type gaussian_impl_128;
    gaussian_impl_type gaussian_impl_256;
    gaussian_impl_type gaussian_impl_512;
    cuda::function<void (float4*, uint, uint, dsfloat, coalesced_vector_type const*, dsfloat const*, dsfloat const*, uint)> shift_rescale;
    static boltzmann_wrapper const kernel;
};

} // namespace mdsim
} // namespace gpu
} // namespace velocities
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_VELOCITIES_BOLTZMANN_KERNEL_HPP */
