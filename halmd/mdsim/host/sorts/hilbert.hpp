/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_SORT_HILBERT_HPP
#define HALMD_MDSIM_HOST_SORT_HILBERT_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/binning.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace sorts {

template <int dimension, typename float_type>
class hilbert
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::binning<dimension, float_type> binning_type;

    static void luaopen(lua_State* L);

    hilbert(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<binning_type> binning
      , std::shared_ptr<logger> logger = std::make_shared<logger>()
    );
    void order();

    connection on_order(std::function<void ()> const& slot)
    {
        return on_order_.connect(slot);
    }

private:
    typedef typename binning_type::cell_size_type cell_size_type;
    typedef typename binning_type::cell_list cell_list;
    typedef typename binning_type::array_type cell_array_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        utility::profiler::accumulator_type order;
        utility::profiler::accumulator_type map;
    };

    unsigned int map(vector_type r, unsigned int depth);

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<box_type const> box_;
    std::shared_ptr<binning_type> binning_;

    /** 1-dimensional Hilbert curve mapping of cell lists */
    std::vector<cell_list const*> map_;
    /** signal emitted after particle ordering */
    signal<void ()> on_order_;
    /** module logger */
    std::shared_ptr<logger> logger_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace sorts
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_SORT_HILBERT_HPP */
