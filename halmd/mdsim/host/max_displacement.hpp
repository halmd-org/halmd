/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_MAX_DISPLACEMENT_HPP
#define HALMD_MDSIM_HOST_MAX_DISPLACEMENT_HPP

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <lua.hpp>
#include <memory>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class max_displacement
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;


    static void luaopen(lua_State* L);

    max_displacement(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<box_type const> box
    );
    void zero();
    float_type compute();

private:
    typedef typename particle_type::size_type size_type;
    typedef typename particle_type::position_array_type position_array_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type zero;
        accumulator_type compute;
    };

    //! system state
    std::shared_ptr<particle_type const> particle_;
    //! simulation box
    std::shared_ptr<box_type const> box_;
    /* particle positions at last neighbour list update */
    std::vector<vector_type> r0_;
    /** cache observer of position updates */
    cache<> position_cache_;
    /** the last calculated displacement */
    float_type displacement_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_MAX_DISPLACEMENT_HPP */
