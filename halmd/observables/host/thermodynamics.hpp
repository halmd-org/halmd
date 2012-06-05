/*
 * Copyright © 2010-2012  Felix Höfling and Peter Colberg
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

#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/data_cache.hpp>
#include <halmd/utility/profiler.hpp>

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
    typedef mdsim::clock clock_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_group particle_group_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    thermodynamics(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<clock_type const> clock
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

    /**
     * Clear cache.
     */
    virtual void clear_cache();

private:
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

    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** particle group */
    std::shared_ptr<particle_group_type> group_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** cached results */
    data_cache<double> en_kin_;
    data_cache<vector_type> v_cm_;
    data_cache<double> en_pot_;
    data_cache<double> virial_;
    data_cache<double> hypervirial_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
