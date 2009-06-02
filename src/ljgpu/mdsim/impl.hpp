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

#ifndef LJGPU_MDSIM_IMPL_HPP
#define LJGPU_MDSIM_IMPL_HPP

#include <boost/type_traits/integral_constant.hpp>

namespace ljgpu
{

struct mdsim_impl_base
{
    // set all implementation properties to false by default
    typedef boost::false_type impl_energy_gpu_sample;
    typedef boost::false_type impl_fixed_size_cell_lists;
    typedef boost::false_type impl_gpu;
    typedef boost::false_type impl_hardsphere_cell_lists;
    typedef boost::false_type impl_hardsphere_event_lists;
    typedef boost::false_type impl_hardsphere_potential;
    typedef boost::false_type impl_host;
    typedef boost::false_type impl_lennard_jones_potential;
    typedef boost::false_type impl_neighbour_lists;
    typedef boost::false_type impl_thermostat;
    typedef boost::false_type impl_trajectory_gpu_sample;
};

struct ljfluid_impl_base : mdsim_impl_base
{
    typedef boost::true_type impl_lennard_jones_potential;
    typedef boost::true_type impl_thermostat;
};

struct ljfluid_impl_gpu_base : mdsim_impl_base
{
    typedef boost::true_type impl_gpu;
    typedef boost::true_type impl_lennard_jones_potential;
    typedef boost::true_type impl_thermostat;
};

struct ljfluid_impl_gpu_square : mdsim_impl_base
{
    typedef boost::true_type impl_energy_gpu_sample;
    typedef boost::true_type impl_gpu;
    typedef boost::true_type impl_lennard_jones_potential;
    typedef boost::true_type impl_thermostat;
    typedef boost::true_type impl_trajectory_gpu_sample;
};

struct ljfluid_impl_gpu_cell : mdsim_impl_base
{
    typedef boost::true_type impl_energy_gpu_sample;
    typedef boost::true_type impl_fixed_size_cell_lists;
    typedef boost::true_type impl_gpu;
    typedef boost::true_type impl_lennard_jones_potential;
    typedef boost::true_type impl_thermostat;
};

struct ljfluid_impl_gpu_neighbour : mdsim_impl_base
{
    typedef boost::true_type impl_energy_gpu_sample;
    typedef boost::true_type impl_fixed_size_cell_lists;
    typedef boost::true_type impl_gpu;
    typedef boost::true_type impl_lennard_jones_potential;
    typedef boost::true_type impl_neighbour_lists;
    typedef boost::true_type impl_thermostat;
    typedef boost::true_type impl_trajectory_gpu_sample;
};

struct ljfluid_impl_host : mdsim_impl_base
{
    typedef boost::true_type impl_host;
    typedef boost::true_type impl_lennard_jones_potential;
    typedef boost::true_type impl_neighbour_lists;
    typedef boost::true_type impl_thermostat;
};

struct hardsphere_impl : mdsim_impl_base
{
    typedef boost::true_type impl_hardsphere_cell_lists;
    typedef boost::true_type impl_hardsphere_event_lists;
    typedef boost::true_type impl_hardsphere_potential;
    typedef boost::true_type impl_host;
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_IMPL_HPP */
