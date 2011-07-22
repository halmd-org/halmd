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

#include <boost/make_shared.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/binning.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace sorts {

template <int dimension, typename float_type>
class hilbert
  : public mdsim::sort<dimension>
{
public:
    typedef mdsim::sort<dimension> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::binning<dimension, float_type> binning_type;
    typedef typename _Base::signal_type signal_type;
    typedef typename _Base::slot_function_type slot_function_type;
    typedef logger logger_type;

    static char const* module_name() { return "hilbert"; }

    static void luaopen(lua_State* L);

    hilbert(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<binning_type> binning
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void order();

    virtual connection on_order(slot_function_type const& slot)
    {
        return on_order_.connect(slot);
    }

private:
    typedef typename binning_type::cell_size_type cell_size_type;
    typedef typename binning_type::cell_list cell_list;
    typedef typename binning_type::cell_lists cell_lists;
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type order;
        accumulator_type map;
        accumulator_type permute;
    };

    unsigned int map(vector_type r, unsigned int depth);

    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<binning_type> binning_;

    /** 1-dimensional Hilbert curve mapping of cell lists */
    std::vector<cell_list const*> map_;
    /** signal emitted after particle ordering */
    signal_type on_order_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace sorts
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_SORT_HILBERT_HPP */
