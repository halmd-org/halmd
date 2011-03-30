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

#include <halmd/modules.hpp>

/**
 * Register Lua bindings for HALMD C++ modules
 *
 * Base classes must be registered *before* derived classes, otherwise
 * Luabind will throw an assertion error or cause a segmentation fault.
 */
HALMD_LUA_API int luaopen_libhalmd(lua_State* L)
{
    luaopen_libhalmd_any_converter(L);
    luaopen_libhalmd_h5(L);
    luaopen_libhalmd_logger(L);
    luaopen_libhalmd_options_parser(L);
    luaopen_libhalmd_po(L);
    luaopen_libhalmd_runner(L);
    luaopen_libhalmd_signal(L);
    luaopen_libhalmd_ublas(L);

    luaopen_libhalmd_io_profiling_writer(L);
    luaopen_libhalmd_io_statevars_writer(L);
    luaopen_libhalmd_io_trajectory_reader(L);
    luaopen_libhalmd_io_trajectory_writer(L);

    luaopen_libhalmd_io_profiling_writers_h5md(L);
    luaopen_libhalmd_io_profiling_writers_log(L);
    luaopen_libhalmd_io_statevars_writers_hdf5(L);
    luaopen_libhalmd_io_trajectory_readers_h5md(L);
    luaopen_libhalmd_io_trajectory_readers_halmd_0_1_x(L);
    luaopen_libhalmd_io_trajectory_writers_h5md(L);

    luaopen_libhalmd_mdsim_box(L);
    luaopen_libhalmd_mdsim_core(L);
    luaopen_libhalmd_mdsim_force(L);
    luaopen_libhalmd_mdsim_integrator(L);
    luaopen_libhalmd_mdsim_neighbour(L);
    luaopen_libhalmd_mdsim_particle(L);
    luaopen_libhalmd_mdsim_position(L);
    luaopen_libhalmd_mdsim_sort(L);
    luaopen_libhalmd_mdsim_velocity(L);

    luaopen_libhalmd_mdsim_integrators_nvt(L);

    luaopen_libhalmd_mdsim_host_force(L);
    luaopen_libhalmd_mdsim_host_neighbour(L);
    luaopen_libhalmd_mdsim_host_particle(L);
    luaopen_libhalmd_mdsim_host_velocity(L);

    luaopen_libhalmd_mdsim_host_forces_lennard_jones(L);
    luaopen_libhalmd_mdsim_host_forces_morse(L);
    luaopen_libhalmd_mdsim_host_forces_power_law(L);
    luaopen_libhalmd_mdsim_host_forces_smooth(L);
    luaopen_libhalmd_mdsim_host_forces_zero(L);
    luaopen_libhalmd_mdsim_host_integrators_verlet(L);
    luaopen_libhalmd_mdsim_host_integrators_verlet_nvt_andersen(L);
    luaopen_libhalmd_mdsim_host_positions_lattice(L);
    luaopen_libhalmd_mdsim_host_positions_phase_space(L);
    luaopen_libhalmd_mdsim_host_sorts_hilbert(L);
    luaopen_libhalmd_mdsim_host_velocities_boltzmann(L);
    luaopen_libhalmd_mdsim_host_velocities_phase_space(L);

    luaopen_libhalmd_observables_density_mode(L);
    luaopen_libhalmd_observables_phase_space(L);
    luaopen_libhalmd_observables_sampler(L);
    luaopen_libhalmd_observables_samples_density_mode(L);
    luaopen_libhalmd_observables_ssf(L);
    luaopen_libhalmd_observables_thermodynamics(L);
    luaopen_libhalmd_observables_utility_wavevector(L);

    luaopen_libhalmd_observables_host_density_mode(L);
    luaopen_libhalmd_observables_host_mean_quartic_displacement(L);
    luaopen_libhalmd_observables_host_mean_square_displacement(L);
    luaopen_libhalmd_observables_host_phase_space(L);
    luaopen_libhalmd_observables_host_samples_phase_space(L);
    luaopen_libhalmd_observables_host_thermodynamics(L);
    luaopen_libhalmd_observables_host_velocity_autocorrelation(L);

    luaopen_libhalmd_random_random(L);
    luaopen_libhalmd_random_host_random(L);

    luaopen_libhalmd_utility_profiler(L);

#ifdef WITH_CUDA
    luaopen_libhalmd_mdsim_gpu_force(L);
    luaopen_libhalmd_mdsim_gpu_neighbour(L);
    luaopen_libhalmd_mdsim_gpu_particle(L);
    luaopen_libhalmd_mdsim_gpu_velocity(L);

    luaopen_libhalmd_mdsim_gpu_forces_lennard_jones(L);
    luaopen_libhalmd_mdsim_gpu_forces_morse(L);
    luaopen_libhalmd_mdsim_gpu_forces_zero(L);
    luaopen_libhalmd_mdsim_gpu_integrators_verlet(L);
    luaopen_libhalmd_mdsim_gpu_integrators_verlet_nvt_andersen(L);
    luaopen_libhalmd_mdsim_gpu_positions_lattice(L);
    luaopen_libhalmd_mdsim_gpu_positions_phase_space(L);
    luaopen_libhalmd_mdsim_gpu_sorts_hilbert(L);
    luaopen_libhalmd_mdsim_gpu_velocities_boltzmann(L);
    luaopen_libhalmd_mdsim_gpu_velocities_phase_space(L);

    luaopen_libhalmd_observables_gpu_density_mode(L);
    luaopen_libhalmd_observables_gpu_phase_space(L);
    luaopen_libhalmd_observables_gpu_samples_phase_space(L);
    luaopen_libhalmd_observables_gpu_thermodynamics(L);

    luaopen_libhalmd_random_gpu_random(L);

    luaopen_libhalmd_utility_gpu_device(L);
#endif /* WITH_CUDA */

    return 0;
}
