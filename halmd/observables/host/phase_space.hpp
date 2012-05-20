/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_HOST_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_HOST_PHASE_SPACE_HPP

#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/particle_group.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
class phase_space
{
public:
    typedef samples::phase_space<dimension, float_type> sample_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::particle_group<particle_type> particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;

    /**
     * Construct phase_space sampler from particle group.
     */
    phase_space(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<particle_group_type const> particle_group
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /**
     * Acquire phase_space sample.
     */
    boost::shared_ptr<sample_type const> acquire();

    /**
     * Set particles from phase_space sample.
     */
    void set(boost::shared_ptr<sample_type const> sample);

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** particle instance to particle group */
    boost::shared_ptr<particle_type> particle_;
    /** particle group */
    boost::shared_ptr<particle_group_type const> particle_group_;
    /** simulation box */
    boost::shared_ptr<box_type const> box_;
    /** simulation clock */
    boost::shared_ptr<clock_type const> clock_;
    /** logger instance */
    boost::shared_ptr<logger_type> logger_;
    /** cached phase_space sample */
    boost::shared_ptr<sample_type> sample_;

    typedef typename sample_type::vector_type vector_type;
    typedef halmd::utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
        accumulator_type reset;
        accumulator_type set;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace host
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_PHASE_SPACE_HPP */
