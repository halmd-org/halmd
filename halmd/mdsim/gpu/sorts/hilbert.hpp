/*
 * Copyright © 2008-2010  Peter Colberg
 * Copyright © 2013       Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_SORTS_HILBERT_HPP
#define HALMD_MDSIM_GPU_SORTS_HILBERT_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/sorts/hilbert_kernel.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace sorts {

template <int dimension, typename float_type>
class hilbert
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename particle_type::gpu_vector_type gpu_vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef hilbert_wrapper<dimension> wrapper_type;

    static void luaopen(lua_State* L);

    hilbert(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );
    void order();

    connection on_order(std::function<void ()> const& slot)
    {
        return on_order_.connect(slot);
    }

private:
    typedef typename particle_type::position_array_type position_array_type;

    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        utility::profiler::accumulator_type order;
        utility::profiler::accumulator_type map;
    };

    void map(cuda::memory::device::vector<unsigned int>& g_map);
    void permutation(cuda::memory::device::vector<unsigned int>& g_map, cuda::memory::device::vector<unsigned int>& g_index);

    std::shared_ptr<particle_type> particle_;
    /** simulation box */
    std::shared_ptr<box_type const> box_;
    /** recursion depth */
    unsigned int depth_;
    /** signal emitted after particle ordering */
    signal<void ()> on_order_;
    /** module logger */
    std::shared_ptr<logger> logger_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace sorts
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_SORTS_HILBERT_HPP */
