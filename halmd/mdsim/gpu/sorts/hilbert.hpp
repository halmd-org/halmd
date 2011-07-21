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

#ifndef HALMD_MDSIM_GPU_SORTS_HILBERT_HPP
#define HALMD_MDSIM_GPU_SORTS_HILBERT_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/sorts/hilbert_kernel.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace sorts {

template <int dimension, typename float_type>
class hilbert
  : public mdsim::sort<dimension>
{
public:
    typedef mdsim::sort<dimension> _Base;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename particle_type::gpu_vector_type gpu_vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef hilbert_wrapper<dimension> wrapper_type;
    typedef typename _Base::signal_type signal_type;
    typedef typename _Base::slot_function_type slot_function_type;
    typedef typename _Base::connection_type connection_type;
    typedef logger logger_type;

    static char const* module_name() { return "hilbert"; }

    static void luaopen(lua_State* L);

    hilbert(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void order();

    virtual connection_type on_order(slot_function_type const& slot)
    {
        return on_order_.connect(slot);
    }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type map;
        accumulator_type permutation;
        accumulator_type order;
    };

    void map(cuda::vector<unsigned int>& g_map);
    void permutation(cuda::vector<unsigned int>& g_map, cuda::vector<unsigned int>& g_index);
    void order(cuda::vector<unsigned int> const& g_index);

    boost::shared_ptr<particle_type> particle_;
    /** recursion depth */
    unsigned int depth_;
    /** signal emitted after particle ordering */
    signal_type on_order_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace gpu
} // namespace sorts
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_SORTS_HILBERT_HPP */
