/*
 * Copyright © 2017  Daniel Kirchner
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

#ifndef HALMD_UTILITY_GPU_TEXTURE_HPP
#define HALMD_UTILITY_GPU_TEXTURE_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/blas/blas.hpp>
#ifdef __CUDACC__
#include <halmd/utility/gpu/texture.cuh>
#else
#include <type_traits>
#endif

namespace cuda {
namespace halmd {

using namespace ::halmd;

template<typename T>
class texture : public cuda::texture<T> {
public:
#ifdef __CUDACC__
    texture(::texture<T> &tex) : cuda::texture<T>(tex) {}
#else /* ! __CUDACC__ */
    void bind(cuda::vector<T> const& array) const {
        cuda::texture<T>::bind(array);
    }
    template<typename U>
    void bind(dsfloat_vector<U> const& array, typename std::enable_if<(std::is_same<U, float>::value || std::is_same<U, float2>::value || std::is_same<U, float4>::value) && std::is_same<U, T>::value>::type* = 0) const {
        cuda::texture<T>::bind(array.storage());
    }
#endif /* ! __CUDACC__ */
};

template<size_t N>
class texture<fixed_vector<dsfloat, N> > {
public:

#ifdef __CUDACC__
    texture(::texture<float4>& tex) : tex_(tex) {}
#else /* ! __CUDACC__ */

    void bind(dsfloat_vector<float4> const& array) const {
        tex_.bind(array.storage());
    }

#endif
private:
    cuda::texture<float4> tex_;
};

} // namespace halmd
} // namespace cuda

#endif /* ! HALMD_UTILTIY_GPU_TEXTURE_HPP */