/*
 * Copyright © 2017       Jake Atwell
 * Copyright © 2016       Manuel Dibak
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_ORIENTATIONS_UNIFORM_HPP
#define HALMD_MDSIM_HOST_ORIENTATIONS_UNIFORM_HPP

#include <lua.hpp>
#include <memory>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace orientations {

template <int dimension, typename float_type>
class uniform
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef random::host::random random_type;

    uniform(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<random_type> random
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
      );

    static void luaopen(lua_State* L);

    void set();

private:
    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type set;
    };

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<random_type> random_;
    std::shared_ptr<logger> logger_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace host
} // namespace orientations
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_ORIENTATIONS_UNIFORM_HPP */
