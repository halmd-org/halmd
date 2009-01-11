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

#include <cuda/vector_types.h>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/sample/sample.hpp>

namespace ljgpu
{

template <typename mdsim_impl>
struct mdsim_traits;

template <>
struct mdsim_traits<ljfluid_impl_gpu_base<2> >
{
    typedef float float_type;
    typedef vector<float, 2> vector_type;
    typedef float2 gpu_vector_type;
    typedef trajectory_gpu_sample<vector_type> sample_type;
    enum { dimension = 2 };
};

template <>
struct mdsim_traits<ljfluid_impl_gpu_base<3> >
{
    typedef float float_type;
    typedef vector<float, 3> vector_type;
    typedef float4 gpu_vector_type;
    typedef trajectory_gpu_sample<vector_type> sample_type;
    enum { dimension = 3 };
};

template <>
struct mdsim_traits<ljfluid_impl_host<2> >
{
    typedef double float_type;
    typedef vector<double, 2> vector_type;
    typedef trajectory_host_sample<vector_type> sample_type;
    enum { dimension = 2 };
};

template <>
struct mdsim_traits<ljfluid_impl_host<3> >
{
    typedef double float_type;
    typedef vector<double, 3> vector_type;
    typedef trajectory_host_sample<vector_type> sample_type;
    enum { dimension = 3 };
};

template <int dimension>
struct mdsim_traits<ljfluid_impl_gpu_square<dimension> >
: public mdsim_traits<ljfluid_impl_gpu_base<dimension> > {};

template <int dimension>
struct mdsim_traits<ljfluid_impl_gpu_cell<dimension> >
: public mdsim_traits<ljfluid_impl_gpu_base<dimension> > {};

template <int dimension>
struct mdsim_traits<ljfluid_impl_gpu_neighbour<dimension> >
: public mdsim_traits<ljfluid_impl_gpu_base<dimension> > {};

template <int dimension>
struct mdsim_traits<hardsphere_impl<dimension> >
: public mdsim_traits<ljfluid_impl_host<dimension> > {};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_TRAITS_HPP */
