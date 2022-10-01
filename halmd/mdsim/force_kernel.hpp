/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2010 Peter Colberg
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

#ifndef HALMD_MDSIM_FORCE_KERNEL_HPP
#define HALMD_MDSIM_FORCE_KERNEL_HPP

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {

/**
 * Diagonal and off-diagonal elements of distance tensor
 */

template <typename float_type>
HALMD_GPU_ENABLED typename type_traits<3, float_type>::stress_tensor_type
make_stress_tensor(fixed_vector<float_type, 3> const& r, fixed_vector<float_type, 3> const& f)
{
    typename type_traits<3, float_type>::stress_tensor_type v;
    v[0] = r[0] * f[0];
    v[1] = r[1] * f[1];
    v[2] = r[2] * f[2];
    v[3] = r[0] * f[1];
    v[4] = r[0] * f[2];
    v[5] = r[1] * f[2];
    return v;
}

template <typename float_type>
HALMD_GPU_ENABLED typename type_traits<3, float_type>::stress_tensor_type
make_stress_tensor(fixed_vector<float_type, 3> const& r)
{
    return make_stress_tensor(r, r);
}

template <typename float_type>
HALMD_GPU_ENABLED typename type_traits<2, float_type>::stress_tensor_type
make_stress_tensor(fixed_vector<float_type, 2> const& r, fixed_vector<float_type, 2> const& f)
{
    typename type_traits<2, float_type>::stress_tensor_type v;
    v[0] = r[0] * f[0];
    v[1] = r[1] * f[1];
    v[2] = r[0] * f[1];
    return v;
}

template <typename float_type>
HALMD_GPU_ENABLED typename type_traits<2, float_type>::stress_tensor_type
make_stress_tensor(fixed_vector<float_type, 2> const& r)
{
    return make_stress_tensor(r, r);
}

/**
 * In GPU memory, the stress tensor contribution from each particle is stored
 * as an array of float_type elements ordered differently than in host memory:
 * row-major order on the host, column-major order on the GPU in the sense that
 * the whole array is considered a (particle number) by
 * (stress_tensor_type::static_size) matrix. [The tensor is actually flattened
 * to a 1-dimensional vector.]
 */

/**
 * write stress tensor in column-major order
 */
template <typename stress_tensor_type>
HALMD_GPU_ENABLED void
write_stress_tensor(float* g_stress, stress_tensor_type const& v, unsigned int stride)
{
    enum { size = stress_tensor_type::static_size };
    for (int i = 0; i < size; ++i) {
        g_stress[i * stride] = v[i];
    }
}

/**
 * read stress tensor in column-major order
 */
template <typename stress_tensor_type>
HALMD_GPU_ENABLED stress_tensor_type
read_stress_tensor(float const* g_stress, unsigned int stride)
{
    stress_tensor_type v;
    enum { size = stress_tensor_type::static_size };
    for (int i = 0; i < size; ++i) {
        v[i] = g_stress[i * stride];
    }
    return v;
}

/**
 * read diagonal elements of stress tensor in column-major order
 */
template <typename vector_type>
HALMD_GPU_ENABLED vector_type
read_stress_tensor_diagonal(float const* g_stress, unsigned int stride)
{
    vector_type v;
    enum { dimension = vector_type::static_size };
    for (int i = 0; i < dimension; ++i) {
        v[i] = g_stress[i * stride];
    }
    return v;
}

#ifdef __CUDACC__

/**
 * read stress tensor in column-major order from texture
 */
template <typename stress_tensor_type>
HALMD_GPU_ENABLED stress_tensor_type
read_stress_tensor(cudaTextureObject_t t_stress_pot, unsigned int i, unsigned int stride)
{
    stress_tensor_type v;
    enum { size = stress_tensor_type::static_size };
    for (int j = 0; j < size; ++j) {
        v[j] = tex1Dfetch<float>(t_stress_pot, i + j * stride);
    }
    return v;
}

/**
 * read diagonal elements of stress tensor in column-major order from texture
 */
template <typename vector_type>
HALMD_GPU_ENABLED vector_type
read_stress_tensor_diagonal(cudaTextureObject_t t_stress_pot, unsigned int i, unsigned int stride)
{
    // as the first d(=dimension) elements are the diagonal elements,
    // read_stress_tensor() can be used for reading the diagonal
    return read_stress_tensor<vector_type>(t_stress_pot, i, stride);
}

#endif // __CUDACC__

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCE_KERNEL_HPP */
