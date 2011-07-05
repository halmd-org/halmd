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

#include <lua.hpp>

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
    typedef utility::profiler profiler_type;
    typedef hilbert_wrapper<dimension> wrapper_type;
    typedef typename _Base::signal_type signal_type;
    typedef typename _Base::slot_function_type slot_function_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type map;
        accumulator_type permutation;
        accumulator_type order;
    };

    static char const* module_name() { return "hilbert"; }

    boost::shared_ptr<particle_type> particle;

    static void luaopen(lua_State* L);

    hilbert(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    void register_runtimes(profiler_type& profiler);
    virtual void order();

    virtual void on_order(slot_function_type const& slot)
    {
        on_order_.connect(slot);
    }

private:
    void map(cuda::vector<unsigned int>& g_map);
    void permutation(cuda::vector<unsigned int>& g_map, cuda::vector<unsigned int>& g_index);
    void order(cuda::vector<unsigned int> const& g_index);

    /** recursion depth */
    unsigned int depth_;
    /** profiling runtime accumulators */
    runtime runtime_;
    /** signal emitted after particle ordering */
    signal_type on_order_;
};

} // namespace mdsim
} // namespace gpu
} // namespace sorts
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_SORTS_HILBERT_HPP */
