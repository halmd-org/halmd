/*
 * Copyright © 2025 Giorgia Marcelli
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

//  rescale the system to a target total energy at specific intervals

#ifndef HALMD_MDSIM_GPU_VELOCITIES_RESCALE_HPP
#define HALMD_MDSIM_GPU_VELOCITIES_RESCALE_HPP

#include <lua.hpp>
#include <memory>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/velocities/rescale_kernel.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/utility/gpu/dsfloat_cuda_vector.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/observables/gpu/thermodynamics.hpp> // for energy reading

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <int dimension, typename float_type>
class rescale
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef observables::gpu::thermodynamics<dimension, float_type> thermo_type;

    // constructor 
    rescale(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<thermo_type> thermo  // use to fetch kinetic/potential energy
      , double target_energy
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    // initialize
    void set();

    // return energy targeted and rate
    void set_target_energy(double energy) { target_energy_ = energy; }
    double target_energy() const { return target_energy_; }

    // bind class to lua
    static void luaopen(lua_State* L);

private:
    // I do not know if this is needed. 
//    typedef typename particle_type::velocity_array_type velocity_array_type;
//    typedef typename particle_type::vector_type vector_type;
//    typedef typename particle_type::gpu_vector_type gpu_vector_type;
    typedef rescale_wrapper<dimension, float_type> wrapper_type;

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<thermo_type> thermo_;

    float_type target_energy_;

    std::shared_ptr<logger> logger_;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime {
        accumulator_type set;
    };

    runtime runtime_;
};

} // namespace velocities
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif // HALMD_MDSIM_GPU_VELOCITIES_RESCALE_HPP
