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

#ifndef LJGPU_MDSIM_TRAITS_HPP
#define LJGPU_MDSIM_TRAITS_HPP

#ifdef WITH_CUDA
# include <vector_types.h>
#endif
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/math/vector4d.hpp>

namespace ljgpu
{

template <typename mdsim_impl, int dimension_>
struct mdsim_traits;

template <int dimension_>
struct mdsim_traits<mdsim_impl_base, dimension_>
{
    enum { dimension = dimension_ };
    typedef energy_sample<dimension> energy_sample_type;
};

#ifdef WITH_CUDA
template <>
struct mdsim_traits<ljfluid_impl_gpu_base, 2> : mdsim_traits<mdsim_impl_base, 2>
{
    typedef float float_type;
    typedef vector<float_type, 2> vector_type;
    typedef float2 gpu_vector_type;
    typedef std::vector<trajectory_host_sample<float_type, 2> > host_sample_type;
    typedef std::vector<trajectory_gpu_sample<2> > gpu_sample_type;
};

template <>
struct mdsim_traits<ljfluid_impl_gpu_base, 3> : mdsim_traits<mdsim_impl_base, 3>
{
    typedef float float_type;
    typedef vector<float_type, 3> vector_type;
    typedef float4 gpu_vector_type;
    typedef std::vector<trajectory_host_sample<float_type, 3> > host_sample_type;
    typedef std::vector<trajectory_gpu_sample<3> > gpu_sample_type;
};

template <int dimension>
struct mdsim_traits<ljfluid_impl_gpu_square, dimension> : mdsim_traits<ljfluid_impl_gpu_base, dimension> {};

template <int dimension>
struct mdsim_traits<ljfluid_impl_gpu_cell, dimension> : mdsim_traits<ljfluid_impl_gpu_base, dimension> {};

template <int dimension>
struct mdsim_traits<ljfluid_impl_gpu_neighbour, dimension> : mdsim_traits<ljfluid_impl_gpu_base, dimension> {};
#endif /* WITH_CUDA */

template <int dimension>
struct mdsim_traits<ljfluid_impl_host, dimension> : mdsim_traits<mdsim_impl_base, dimension>
{
    typedef double float_type;
    typedef vector<float_type, dimension> vector_type;
    typedef std::vector<trajectory_host_sample<float_type, dimension> > host_sample_type;
};

template <int dimension>
struct mdsim_traits<hardsphere_impl, dimension> : mdsim_traits<ljfluid_impl_host, dimension> {};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_TRAITS_HPP */
