/*
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

#ifndef HALMD_MDSIM_GPU_ORIENTATIONS_UNIFORM_HPP
#define HALMD_MDSIM_GPU_ORIENTATIONS_UNIFORM_HPP

#include <lua.hpp>
#include <memory>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {

template <int dimension, typename float_type, typename RandomNumberGenerator>
class uniform
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef random::gpu::random<RandomNumberGenerator> random_type;
    typedef typename RandomNumberGenerator::rng_type rng_type;
    typedef typename type_traits<dimension, float>::vector_type gpu_vector_type;
    typedef typename type_traits<dimension, unsigned int>::vector_type index_type;

    uniform(std::shared_ptr<particle_type> particle, std::shared_ptr<random_type> random);// : particle_(particle) { }
    static void luaopen(lua_State* L);

    void set();

private:
    //typedef typename particle_type::position_array_type position_array_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type set;
    };

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<random_type> random_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_ORIENTATIONS_UNIFORM_HPP */
