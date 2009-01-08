/* Lennard-Jones fluid simulation using CUDA
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef LJGPU_LJFLUID_LJFLUID_TRAITS_HPP
#define LJGPU_LJFLUID_LJFLUID_TRAITS_HPP

#include <cuda/vector_types.h>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/sample/sample.hpp>

namespace ljgpu
{

template <int dimension>
struct ljfluid_gpu_traits;

template <>
struct ljfluid_gpu_traits<2>
{
    typedef float float_type;
    typedef vector<float, 2> vector_type;
    typedef float2 gpu_vector_type;
    typedef trajectory_gpu_sample<vector_type> trajectory_sample;
};

template <>
struct ljfluid_gpu_traits<3>
{
    typedef float float_type;
    typedef vector<float, 3> vector_type;
    typedef float4 gpu_vector_type;
    typedef trajectory_gpu_sample<vector_type> trajectory_sample;
};

template <int dimension>
struct ljfluid_host_traits;

template <>
struct ljfluid_host_traits<2>
{
    typedef double float_type;
    typedef vector<double, 2> vector_type;
    typedef trajectory_host_sample<vector_type> trajectory_sample;
};

template <>
struct ljfluid_host_traits<3>
{
    typedef double float_type;
    typedef vector<double, 3> vector_type;
    typedef trajectory_host_sample<vector_type> trajectory_sample;
};

} // namespace ljgpu

#endif /* ! LJGPU_LJFLUID_LJFLUID_TRAITS_HPP */
