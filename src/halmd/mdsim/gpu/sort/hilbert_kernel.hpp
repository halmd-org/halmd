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

#ifndef HALMD_MDSIM_GPU_HILBERT_KERNEL_HPP
#define HALMD_MDSIM_GPU_HILBERT_KERNEL_HPP

#include <cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension>
struct hilbert_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;

    /** Hilbert space-filling curve recursion depth */
    cuda::symbol<float> depth;
    /** cubic box edgle length */
    cuda::symbol<vector_type> box_length;
    /** generate Hilbert space-filling curve */
    cuda::function<void (float4 const*, unsigned int*)> map;

    static hilbert_wrapper const kernel;
};

template <int dimension>
hilbert_wrapper<dimension> const& get_hilbert_kernel()
{
    return hilbert_wrapper<dimension>::kernel;
}

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_HILBERT_KERNEL_HPP */
