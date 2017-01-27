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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/blas/blas.hpp>
#ifdef __CUDACC__
#include <halmd/utility/gpu/texture.cuh>
#endif

namespace cuda {
namespace halmd {

using namespace ::halmd;

template<typename T>
class texture : public cuda::texture<T> {
public:
#ifdef __CUDACC__
    texture(::texture<T> &tex) : cuda::texture<T>(tex) {}
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
