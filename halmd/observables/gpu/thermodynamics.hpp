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

#ifndef HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/observables/gpu/thermodynamics_kernel.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/data_cache.hpp>
#include <halmd/utility/profiler.hpp>

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
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::clock clock_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    thermodynamics(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
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
    /** module dependencies */
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<particle_type const> particle_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

    /** cached results */
    data_cache<double> en_kin_;
    data_cache<vector_type> v_cm_;
    data_cache<double> en_pot_;
    data_cache<double> virial_;
    data_cache<double> hypervirial_;

    /** functor for total kinetic energy */
    reduction<gpu::kinetic_energy<dimension, dsfloat> > compute_en_kin_;
    /** functor for centre-of-mass velocity */
    reduction<gpu::velocity_of_centre_of_mass<dimension, dsfloat> > compute_v_cm_;
    /** functor for total potential energy */
    reduction<gpu::potential_energy<dsfloat> > compute_en_pot_;
    /** functor for total virial sum */
    reduction<gpu::virial<dimension, dsfloat> > compute_virial_;

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

} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
