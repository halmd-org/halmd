/*
 * Copyright © 2010-2023 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2010-2012 Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/profiler.hpp>

#include <lua.hpp>

#include <functional>
#include <memory>
#include <tuple>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
class thermodynamics
  : public observables::thermodynamics<dimension>
{
private:
    typedef observables::thermodynamics<dimension> _Base;

public:
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::stress_tensor_type stress_tensor_type;
    typedef mdsim::box<dimension> box_type;
    typedef std::function<double ()> volume_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_group particle_group_type;

    static void luaopen(lua_State* L);

    thermodynamics(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , volume_type volume = nullptr
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Returns number of particles.
     */
    virtual unsigned int particle_number() const;

    /**
     * Returns the volume obtained from calling volume_().
     */
    virtual double volume() const;

    /**
     * Compute mean kinetic energy per particle.
     */
    virtual double en_kin();

    /**
     * Compute total force.
     */
    virtual vector_type const& total_force();

    /**
     * Compute velocity of centre of mass
     */
    virtual vector_type const& v_cm();

    /**
     * Compute centre of mass
     */
    virtual vector_type const& r_cm();

    /**
     * Compute mean particle mass.
     */
    virtual double mean_mass();

    /**
     * Compute mean potential energy per particle.
     */
    virtual double en_pot();

    /**
     * Compute mean virial per particle.
     */
    virtual double virial();


    /**
     * Compute mean stress tensor elements per particle.
     */
    virtual stress_tensor_type const& stress_tensor();

private:
    typedef typename particle_type::size_type size_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::force_array_type force_array_type;
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;
    typedef typename particle_group_type::array_type group_array_type;

    /** system state */
    std::shared_ptr<particle_type> particle_;
    /** particle group */
    std::shared_ptr<particle_group_type> group_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** reference volume */
    volume_type volume_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** mean kinetic energy per particle */
    double en_kin_;
    /** total force */
    vector_type force_;
    /** velocity of centre of mass */
    vector_type v_cm_;
    /** centre of mass */
    vector_type r_cm_;
    /** mean mass */
    double mean_mass_;
    /** mean potential energy per particle */
    double en_pot_;
    /** mean virial per particle */
    double virial_;
    /** mean stress tensor elements per particle */
    stress_tensor_type stress_tensor_;

    /** cache observers of mean kinetic energy per particle */
    std::tuple<cache<>, cache<>> en_kin_cache_;
    /** cache observers of total force */
    std::tuple<cache<>, cache<>> force_cache_;
    /** cache observers of velocity of centre of mass */
    std::tuple<cache<>, cache<>> v_cm_cache_;
    /** cache observers of centre of mass */
    std::tuple<cache<>, cache<>> r_cm_cache_;
    /** cache observers of mean potential energy per particle */
    std::tuple<cache<>, cache<>> en_pot_cache_;
    /** cache observers of mean virial per particle */
    std::tuple<cache<>, cache<>> virial_cache_;
    /** cache observers of mean stress tensor elements per particle */
    std::tuple<cache<>, cache<>, cache<>> stress_tensor_cache_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type en_kin;
        accumulator_type force;
        accumulator_type v_cm;
        accumulator_type r_cm;
        accumulator_type en_pot;
        accumulator_type virial;
        accumulator_type stress_tensor;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
