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

#ifndef HALMD_MDSIM_GPU_POSITIONS_LATTICE_HPP
#define HALMD_MDSIM_GPU_POSITIONS_LATTICE_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {

template <int dimension, typename float_type>
class lattice
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef logger logger_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename type_traits<dimension, float>::vector_type gpu_vector_type;
    typedef typename type_traits<dimension, unsigned int>::vector_type index_type;

    static void luaopen(lua_State* L);

    lattice(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , typename box_type::vector_type const& slab
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    void set();

    typename box_type::vector_type const& slab() const { return slab_; }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type set;
    };

    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<logger_type> logger_;
    /** slab extents for each direction as fraction of the edge length of the box */
    typename box_type::vector_type slab_;

    /**
     *  assign range of particle positions [first, last) to fcc lattice
     *  of extents 'length' with origin at 'offset'
     */
    template <typename position_iterator>
    void fcc(
        position_iterator first, position_iterator last
      , gpu_vector_type const& length, gpu_vector_type const& offset
    );

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITIONS_LATTICE_HPP */
