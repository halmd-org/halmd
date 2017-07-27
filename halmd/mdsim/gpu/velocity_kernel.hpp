/*
 * Copyright © 2008-2010, 2012 Peter Colberg
 * Copyright © 2010 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_VELOCITY_KERNEL_HPP
#define HALMD_MDSIM_GPU_VELOCITY_KERNEL_HPP

#include <boost/mpl/if.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
struct velocity_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename type_traits<4, float_type>::gpu::ptr_type ptr_type;
    cuda::function<void (ptr_type, unsigned int, unsigned int, dsfloat)> rescale;
    cuda::function<void (ptr_type, unsigned int const*, unsigned int, unsigned int, dsfloat)> rescale_group;
    cuda::function<void (ptr_type, unsigned int, unsigned int, fixed_vector<dsfloat, dimension>)> shift;
    cuda::function<void (ptr_type, unsigned int const*, unsigned int, unsigned int, fixed_vector<dsfloat, dimension>)> shift_group;
    cuda::function<void (ptr_type, unsigned int, unsigned int, fixed_vector<dsfloat, dimension>, dsfloat)> shift_rescale;
    cuda::function<void (ptr_type, unsigned int const*, unsigned int, unsigned int, fixed_vector<dsfloat, dimension>, dsfloat)> shift_rescale_group;
    static velocity_wrapper const kernel;
};

template <int dimension, typename float_type>
velocity_wrapper<dimension, float_type> const& get_velocity_kernel()
{
    return velocity_wrapper<dimension, float_type>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_VELOCITY_KERNEL_HPP */
