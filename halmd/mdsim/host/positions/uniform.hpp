/*
 * Copyright © 2016 Arthur Straube
 * Copyright © 2016 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_POSITIONS_UNIFORM_HPP
#define HALMD_MDSIM_HOST_POSITIONS_UNIFORM_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/position.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

template <int dimension, typename float_type>
class uniform
  : public position
{
    typedef mdsim::host::position _Base;

public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef halmd::random::host::random rng_type;

    uniform(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<rng_type> rng
      , vector_type const& slab
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Sample random particle positions uniformly within slab centred around origin.
     */
    void set();

    /**
     * Return slab extents as fractions of box edges.
     */
    vector_type const& slab() const { return slab_; }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** system state */
    std::shared_ptr<particle_type> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** random number generator */
    std::shared_ptr<rng_type> rng_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** slab extents for each direction as fraction of the edge lengths of the box */
    vector_type slab_;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type set;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace host
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POSITIONS_UNIFORM_HPP */
