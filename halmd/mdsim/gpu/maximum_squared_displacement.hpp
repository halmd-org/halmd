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

#ifndef HALMD_MDSIM_GPU_MAXIMUM_SQUARED_DISPLACEMENT_HPP
#define HALMD_MDSIM_GPU_MAXIMUM_SQUARED_DISPLACEMENT_HPP

#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/maximum_squared_displacement_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * Compute maximum squared displacement
 */
template <int dimension, typename float_type>
class maximum_squared_displacement
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;

    static void luaopen(lua_State* L);

    maximum_squared_displacement(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
    );
    void zero();
    float_type compute();

private:
    typedef typename particle_type::vector_type vector_type;
    typedef typename maximum_squared_displacement_wrapper<dimension>::displacement_impl_type displacement_impl_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type zero;
        accumulator_type compute;
    };

    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<box_type const> box_;

    cuda::config dim_reduce_;
    static displacement_impl_type get_displacement_impl(int threads);
    displacement_impl_type displacement_impl_;

    /** particle positions at last maximum displacement list update */
    cuda::vector<float4> g_r0_;
    /** block-reduced squared particle distances */
    cuda::vector<float> g_rr_;
    /** block-reduced squared particle distances */
    cuda::host::vector<float> h_rr_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_MAXIMUM_SQUARED_DISPLACEMENT_HPP */
