/*
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2010-2012 Peter Colberg
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

#ifndef HALMD_OBSERVABLES_HOST_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_HOST_THERMODYNAMICS_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/profiler.hpp>

#include <lua.hpp>

#include <memory>
#include <tuple>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
class thermodynamics
  : public observables::thermodynamics<dimension>
{
private:
    typedef observables::thermodynamics<dimension> _Base;

public:
    typedef typename _Base::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::force<dimension, float_type> force_type;
    typedef mdsim::host::particle_group particle_group_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    thermodynamics(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<force_type> force
      , std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /**
     * Returns number of particles.
     */
    virtual unsigned int nparticle() const;

    /**
     * Returns volume.
     */
    virtual double volume() const;

    /**
     * Compute mean kinetic energy per particle.
     */
    virtual double en_kin();

    /**
     * Compute mean potential energy per particle.
     */
    virtual vector_type const& v_cm();

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
     * Compute mean hypervirial per particle.
     */
    virtual double hypervirial();

private:
    typedef typename particle_type::size_type size_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::velocity_type velocity_type;
    typedef typename particle_type::mass_array_type mass_array_type;
    typedef typename particle_type::mass_type mass_type;
    typedef typename force_type::en_pot_array_type en_pot_array_type;
    typedef typename force_type::stress_pot_array_type stress_pot_array_type;
    typedef typename force_type::hypervirial_array_type hypervirial_array_type;
    typedef typename particle_group_type::array_type group_array_type;

    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** particle forces */
    std::shared_ptr<force_type> force_;
    /** particle group */
    std::shared_ptr<particle_group_type> group_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** mean kinetic energy per particle */
    double en_kin_;
    /** velocity of centre of mass */
    vector_type v_cm_;
    /** mean mass */
    double mean_mass_;
    /** mean potential energy per particle */
    double en_pot_;
    /** mean virial per particle */
    double virial_;
    /** mean hypervirial per particle */
    double hypervirial_;

    /** cache observers of mean kinetic energy per particle */
    std::tuple<cache<>, cache<>, cache<>> en_kin_cache_;
    /** cache observers of mean potential energy per particle */
    std::tuple<cache<>, cache<>, cache<>> v_cm_cache_;
    /** cache observers of mean potential energy per particle */
    std::tuple<cache<>, cache<>> en_pot_cache_;
    /** cache observers of mean virial per particle */
    std::tuple<cache<>, cache<>> virial_cache_;
    /** cache observers of mean hypervirial per particle */
    std::tuple<cache<>, cache<>> hypervirial_cache_;

    typedef halmd::utility::profiler profiler_type;
    typedef profiler_type::accumulator_type accumulator_type;
    typedef profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type en_kin;
        accumulator_type v_cm;
        accumulator_type en_pot;
        accumulator_type virial;
        accumulator_type hypervirial;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
