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

#ifndef HALMD_OBSERVABLES_GPU_PHASE_SPACE_KERNEL_CUH
#define HALMD_OBSERVABLES_GPU_PHASE_SPACE_KERNEL_CUH

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension>
struct phase_space_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    /** positions, types */
    cuda::texture<float4> r;
    /** minimum image vectors */
    cuda::texture<coalesced_vector_type> image;
    /** velocities, tags */
    cuda::texture<float4> v;
    /** sample phase space for all particle of a single species */
    cuda::function<void (unsigned int const*, coalesced_vector_type*, coalesced_vector_type*, unsigned int, vector_type)> sample;

    static phase_space_wrapper const kernel;
};

template <int dimension>
phase_space_wrapper<dimension> const& get_phase_space_kernel()
{
    return phase_space_wrapper<dimension>::kernel;
}

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PHASE_SPACE_KERNEL_CUH */
