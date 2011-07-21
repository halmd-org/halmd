/*
 * Copyright © 2011  Felix Höfling
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

#ifndef HALMD_TEST_UNIT_MDSIM_BOX_KERNEL_HPP
#define HALMD_TEST_UNIT_MDSIM_BOX_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>

template <int dimension, typename float_type>
struct box_kernel_wrapper
{
    typedef typename halmd::mdsim::type_traits<dimension, float_type>::vector_type vector_type;
    typedef typename halmd::mdsim::type_traits<dimension, float_type>::gpu::coalesced_vector_type coalesced_vector_type;

    cuda::function <void (coalesced_vector_type*, coalesced_vector_type*, const vector_type, const unsigned int)> reduce_periodic;

    static box_kernel_wrapper kernel;
};

#endif /* ! HALMD_TEST_UNIT_MDSIM_BOX_KERNEL_HPP */
