/*
 * Copyright © 2011  Peter Colberg
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

#include <cuda.h> // CUDA_VERSION

// boost::mpl::sizeof_ compiles only with CUDA 3.2 or higher
#if CUDA_VERSION >= 3020

#include <boost/mpl/sizeof.hpp>
#include <boost/type_traits/alignment_of.hpp>

#include "variant_kernel.hpp"

namespace variant_kernel
{

/**
 * Returns alignment of type in first element of result array.
 *
 * @param result Result array with one element.
 */
template <typename T>
__global__ void alignment_of(size_t* result)
{
    *result = boost::alignment_of<T>::value;
}

/**
 * Returns size in bytes of type in first element of result array.
 *
 * @param result Result array with one element.
 */
template <typename T>
__global__ void sizeof_(size_t* result)
{
    *result = boost::mpl::sizeof_<T>::value;
}

} // namespace variant_kernel

// wrap test kernels for arbitrary type
template <typename T>
prerequisites_wrapper<T> prerequisites_wrapper<T>::kernel = {
    &variant_kernel::alignment_of<T>
  , &variant_kernel::sizeof_<T>
};

// instantiate test kernels for specific types
template class prerequisites_wrapper<char>;
template class prerequisites_wrapper<short>;
template class prerequisites_wrapper<int>;
template class prerequisites_wrapper<long>;
template class prerequisites_wrapper<float>;
template class prerequisites_wrapper<double>;

template class prerequisites_wrapper<float1>;
template class prerequisites_wrapper<float2>;
template class prerequisites_wrapper<float3>;
template class prerequisites_wrapper<float4>;

template class prerequisites_wrapper<vector4<char> >;
template class prerequisites_wrapper<vector4<short> >;
template class prerequisites_wrapper<vector4<int> >;
template class prerequisites_wrapper<vector4<long> >;
template class prerequisites_wrapper<vector4<float> >;
template class prerequisites_wrapper<vector4<double> >;

template class prerequisites_wrapper<vector4<float1> >;
template class prerequisites_wrapper<vector4<float2> >;
template class prerequisites_wrapper<vector4<float3> >;
template class prerequisites_wrapper<vector4<float4> >;

template class prerequisites_wrapper<vector4_align8<char> >;
template class prerequisites_wrapper<vector4_align8<short> >;
template class prerequisites_wrapper<vector4_align8<int> >;
template class prerequisites_wrapper<vector4_align8<long> >;
template class prerequisites_wrapper<vector4_align8<float> >;
template class prerequisites_wrapper<vector4_align8<double> >;

template class prerequisites_wrapper<vector4_align8<float1> >;
template class prerequisites_wrapper<vector4_align8<float2> >;
template class prerequisites_wrapper<vector4_align8<float3> >;
template class prerequisites_wrapper<vector4_align8<float4> >;

template class prerequisites_wrapper<vector4_align16<char> >;
template class prerequisites_wrapper<vector4_align16<short> >;
template class prerequisites_wrapper<vector4_align16<int> >;
template class prerequisites_wrapper<vector4_align16<long> >;
template class prerequisites_wrapper<vector4_align16<float> >;
template class prerequisites_wrapper<vector4_align16<double> >;

template class prerequisites_wrapper<vector4_align16<float1> >;
template class prerequisites_wrapper<vector4_align16<float2> >;
template class prerequisites_wrapper<vector4_align16<float3> >;
template class prerequisites_wrapper<vector4_align16<float4> >;

#endif /* CUDA_VERSION >= 3020 */
