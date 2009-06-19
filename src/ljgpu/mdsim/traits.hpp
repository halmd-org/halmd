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
#include <boost/mpl/and.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/math/vector4d.hpp>

namespace ljgpu
{

template <typename impl, int dimension_, typename Enable = void>
struct mdsim_traits;

#ifdef WITH_CUDA
template <typename impl, int dimension_>
struct mdsim_traits<impl, dimension_, typename boost::enable_if<typename impl::impl_gpu>::type>
{
    enum { dimension = dimension_ };
    typedef energy_sample<dimension> energy_sample_type;
    typedef float float_type;
    typedef vector<float_type, dimension> vector_type;
    typedef typename boost::mpl::if_c<dimension == 3, float4, float2>::type gpu_vector_type;
    typedef std::vector<trajectory_host_sample<float_type, dimension> > host_sample_type;
    typedef std::vector<trajectory_gpu_sample<dimension> > gpu_sample_type;
};
#endif /* WITH_CUDA */

template <typename impl, int dimension_>
struct mdsim_traits<impl, dimension_, typename boost::enable_if<boost::mpl::and_<typename impl::impl_host, typename impl::impl_lennard_jones_potential> >::type>
{
    enum { dimension = dimension_ };
    typedef energy_sample<dimension> energy_sample_type;
#ifdef USE_HOST_SINGLE_PRECISION
    typedef float float_type;
#else
    typedef double float_type;
#endif
    typedef vector<float_type, dimension> vector_type;
    typedef std::vector<trajectory_host_sample<float_type, dimension> > host_sample_type;
};

template <typename impl, int dimension_>
struct mdsim_traits<impl, dimension_, typename boost::enable_if<boost::mpl::and_<typename impl::impl_host, typename impl::impl_hardsphere_potential> >::type>
{
    enum { dimension = dimension_ };
    typedef energy_sample<dimension> energy_sample_type;
    typedef double float_type;
    typedef vector<float_type, dimension> vector_type;
    typedef std::vector<trajectory_host_sample<float_type, dimension> > host_sample_type;
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_TRAITS_HPP */
