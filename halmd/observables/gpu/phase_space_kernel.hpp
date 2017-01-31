/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_PHASE_SPACE_KERNEL_CUH
#define HALMD_OBSERVABLES_GPU_PHASE_SPACE_KERNEL_CUH

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/gpu/texture.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension>
struct phase_space_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    /** positions, types */
    cuda::halmd::texture<float4> r;
    /** minimum image vectors */
    cuda::texture<coalesced_vector_type> image;
    /** sample position for all particle of a single species */
    cuda::function<void (unsigned int const*, float4*, vector_type, unsigned int)> sample_position;
    /** shift particle positions to range (-L/2, L/2) */
    cuda::function<void (unsigned int const*, float4*, coalesced_vector_type*, vector_type, unsigned int)> reduce_periodic;

    static phase_space_wrapper const kernel;
};

namespace detail {

template<typename T>
struct sample_ptr_type {
    typedef T* ptr_type;
};

template<size_t dimension>
struct sample_ptr_type<fixed_vector<dsfloat, dimension> > {
    typedef typename mdsim::type_traits<dimension, dsfloat>::gpu::ptr_type ptr_type;
};

} // namespace detail

template<typename input_data_type, typename sample_data_type = input_data_type>
struct phase_space_sample_wrapper
{
    cuda::halmd::texture<sample_data_type> input;
    typedef typename detail::sample_ptr_type<input_data_type>::ptr_type ptr_type;
    cuda::function<void (unsigned int const*, sample_data_type*, unsigned int)> sample;
    cuda::function<void (unsigned int const*, ptr_type, unsigned int)> set;

    static phase_space_sample_wrapper const kernel;
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
