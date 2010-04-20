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

#ifndef HALMD_MDSIM_GPU_HILBERT_WRAPPER_CUH
#define HALMD_MDSIM_GPU_HILBERT_WRAPPER_CUH

#include <boost/mpl/if.hpp>

#include <cuda_wrapper.hpp>

namespace halmd { namespace mdsim { namespace gpu
{

template <int dimension>
struct hilbert_wrapper
{
    typedef typename boost::mpl::if_c<dimension == 3, float3, float2>::type vector_type;

    /** Hilbert space-filling curve recursion depth */
    static cuda::symbol<float> depth;
    /** cubic box edgle length */
    static cuda::symbol<vector_type> box_length;
    /** generate Hilbert space-filling curve */
    static cuda::function<void (float4 const*, unsigned int*)> map;
};

}}} // namespace halmd::mdsim::gpu

#endif /* ! HALMD_MDSIM_GPU_HILBERT_WRAPPER_CUH */
