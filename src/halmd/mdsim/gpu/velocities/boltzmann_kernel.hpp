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

#ifndef HALMD_MDSIM_GPU_VELOCITIES_BOLTZMANN_KERNEL_HPP
#define HALMD_MDSIM_GPU_VELOCITIES_BOLTZMANN_KERNEL_HPP

#include <boost/mpl/if.hpp>

#include <cuda_wrapper.hpp>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>

//
// Maxwell-Boltzmann distribution at accurate temperature
//

namespace halmd
{
namespace mdsim { namespace gpu { namespace velocities
{

using numeric::gpu::blas::dsfloat;

template <
    int dimension
  , typename float_type
  , typename rng_type
>
struct boltzmann_wrapper
{
    typedef typename boost::mpl::if_c<dimension == 3, float4, float2>::type coalesced_vector_type;
    typedef cuda::function<void (float4*, uint, uint, float, coalesced_vector_type*, dsfloat*)> gaussian_impl_type;
    gaussian_impl_type gaussian_impl_32;
    gaussian_impl_type gaussian_impl_64;
    gaussian_impl_type gaussian_impl_128;
    gaussian_impl_type gaussian_impl_256;
    gaussian_impl_type gaussian_impl_512;
    cuda::function<void (float4*, uint, uint, dsfloat, coalesced_vector_type const*, dsfloat const*, uint)> shift_rescale;
    cuda::symbol<rng_type> rng;
    static boltzmann_wrapper const kernel;
};

}}} // namespace mdsim::gpu::velocities

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_VELOCITIES_BOLTZMANN_KERNEL_HPP */
