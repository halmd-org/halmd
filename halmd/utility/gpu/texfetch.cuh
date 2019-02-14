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
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {

template<typename float_type>
struct texFetch;

template<>
struct texFetch<float>
{
    template<typename T>
    static __device__ T fetch(texture<T> tex, int i)
    {
        return tex1Dfetch(tex, i);
    }
};

template<>
struct texFetch<dsfloat>
{
    template<typename T>
    static __device__ tuple<T, T> fetch(texture<T> tex, int i)
    {
        return make_tuple(tex1Dfetch(tex, i), tex1Dfetch(tex, i + GTDIM));
    }
};

} // namespace halmd
