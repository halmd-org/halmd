/*
 * Copyright © 2021 Jaslo Ziska
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

#ifndef HALMD_TEST_TOOLS_VECTOR_HPP
#define HALMD_TEST_TOOLS_VECTOR_HPP

#include <vector>

#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
#endif

template <typename T>
std::vector<T> get_host_vector(std::vector<T> const& v)
{
    return v;
}

#ifdef HALMD_WITH_GPU
template <typename T>
cuda::host::vector<T> get_host_vector(cuda::vector<T> const& g_v)
{
    cuda::host::vector<T> v(g_v.size());
    cuda::copy(g_v.begin(), g_v.end(), v.begin());
    return v;
}
#endif

#endif // ! HALMD_TEST_TOOLS_VECTOR_HPP
