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

#ifndef HALMD_MDSIM_GPU_VELOCITY_KERNEL_HPP
#define HALMD_MDSIM_GPU_VELOCITY_KERNEL_HPP

#include <boost/mpl/if.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension>
struct velocity_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    cuda::function<void (float4*, unsigned int, dsfloat)> rescale;
    cuda::function<void (float4*, unsigned int, fixed_vector<dsfloat, dimension>)> shift;
    cuda::function<void (float4*, unsigned int, fixed_vector<dsfloat, dimension>, dsfloat)> shift_rescale;
    cuda::symbol<unsigned int> nbox;
    static velocity_wrapper const kernel;
};

template <int dimension>
velocity_wrapper<dimension> const& get_velocity_kernel()
{
    return velocity_wrapper<dimension>::kernel;
}

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_VELOCITY_KERNEL_HPP */
