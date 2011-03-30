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

#ifndef HALMD_MODULES_HPP
#define HALMD_MODULES_HPP

#include <lua.hpp>

#include <halmd/config.hpp>

HALMD_LUA_API int luaopen_libhalmd_any_converter(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_h5(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_logger(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_options_parser(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_po(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_runner(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_signal(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_ublas(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_io_profiling_writer(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_statevars_writer(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_trajectory_reader(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_trajectory_writer(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_io_profiling_writers_h5md(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_profiling_writers_log(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_statevars_writers_hdf5(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_trajectory_readers_h5md(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_trajectory_readers_halmd_0_1_x(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_io_trajectory_writers_h5md(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_mdsim_box(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_core(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_force(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_integrator(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_neighbour(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_particle(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_position(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_sort(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_velocity(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_mdsim_integrators_nvt(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_force(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_neighbour(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocity(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_lennard_jones(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_morse(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_power_law(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_smooth(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_zero(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet_nvt_andersen(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_lattice(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_sorts_hilbert(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocities_boltzmann(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocities_phase_space(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_observables_density_mode(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_sampler(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_samples_density_mode(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_ssf(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_thermodynamics(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_utility_wavevector(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_observables_host_density_mode(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_host_mean_quartic_displacement(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_host_mean_square_displacement(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_host_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_host_samples_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_host_thermodynamics(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_host_velocity_autocorrelation(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_random_random(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_random_host_random(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_utility_profiler(lua_State* L);

#ifdef WITH_CUDA
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_force(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_neighbour(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_particle(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocity(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_forces_lennard_jones(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_forces_morse(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_forces_zero(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet_nvt_andersen(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_positions_lattice(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_positions_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_sorts_hilbert(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocities_boltzmann(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocities_phase_space(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_density_mode(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_gpu_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_phase_space(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_observables_gpu_thermodynamics(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_random_gpu_random(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd_utility_gpu_device(lua_State* L);
#endif /* WITH_CUDA */

#endif /* ! HALMD_MODULES_HPP */
