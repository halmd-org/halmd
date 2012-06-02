/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_POSITIONS_LATTICE_HPP
#define HALMD_MDSIM_HOST_POSITIONS_LATTICE_HPP

#include <lua.hpp>
#include <memory>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

template <int dimension, typename float_type>
class lattice
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    lattice(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , vector_type const& slab
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );
    void set();

    vector_type const& slab() const { return slab_; }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type set;
    };

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<box_type const> box_;
    std::shared_ptr<logger_type> logger_;
    /** slab extents for each direction as fraction of the edge length of the box */
    vector_type slab_;

    /**
     *  assign range of particle positions [first, last) to fcc lattice
     *  of extents 'length' with origin at 'offset'
     */
    template <typename position_iterator>
    void fcc(
        position_iterator first, position_iterator last
      , vector_type const& length, vector_type const& offset
    );

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace host
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POSITIONS_LATTICE_HPP */
