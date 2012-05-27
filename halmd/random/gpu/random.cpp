/*
 * Copyright © 2010-2011  Peter Colberg and Felix Höfling
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

#include <boost/shared_ptr.hpp>

#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace random {
namespace gpu {

template <typename RandomNumberGenerator>
random<RandomNumberGenerator>::random(
    unsigned int seed
  , unsigned int blocks
  , unsigned int threads
  , unsigned int shuffle_threads
)
  // allocate random number generator state
  : rng_(blocks, threads)
  , shuffle_threads_(shuffle_threads)
{
    random::seed(seed);
}

template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::seed(unsigned int seed)
{
    rng_.seed(seed);
}

/**
 * fill array with uniform random numbers in [0.0, 1.0)
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::uniform(cuda::vector<float>& g_v)
{
    cuda::configure(rng_.dim.grid, rng_.dim.block);
    get_random_kernel<rng_type>().uniform(g_v, g_v.size(), rng_.rng());
    cuda::thread::synchronize();
}

/**
 * fill array with random integers in [0, 2^32-1]
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::get(cuda::vector<unsigned int>& g_v)
{
    cuda::configure(rng_.dim.grid, rng_.dim.block);
    get_random_kernel<rng_type>().get(g_v, g_v.size(), rng_.rng());
    cuda::thread::synchronize();
}

/**
 * fill array with normal distributed random numbers in [0.0, 1.0)
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::normal(cuda::vector<float>& g_v, float mean, float sigma)
{
    cuda::configure(rng_.dim.grid, rng_.dim.block);
    get_random_kernel<rng_type>().normal(g_v, g_v.size(), mean, sigma, rng_.rng());
    cuda::thread::synchronize();
}

template <typename RandomNumberGenerator>
unsigned int random<RandomNumberGenerator>::defaults::blocks() {
    return 32;
}
template <typename RandomNumberGenerator>
unsigned int random<RandomNumberGenerator>::defaults::threads() {
    return 32 << DEVICE_SCALE;
}
template <typename RandomNumberGenerator>
unsigned int random<RandomNumberGenerator>::defaults::shuffle_threads() {
    return 128;
}

template <typename random_type>
static std::function<void (unsigned int)>
wrap_seed(boost::shared_ptr<random_type> self)
{
    return [=](unsigned int seed) {
        self->seed(seed);
    };
}

template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luabind;
    static std::string const class_name = RandomNumberGenerator::name();
    module(L, "libhalmd")
    [
        namespace_("random")
        [
            namespace_("gpu")
            [
                class_<random, boost::shared_ptr<random> >(class_name.c_str())
                    .def(constructor<>())
                    .property("seed", &wrap_seed<random>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_random_gpu_random(lua_State* L)
{
    random<rand48>::luaopen(L);
    return 0;
}

} // namespace random
} // namespace gpu
template class random::gpu::random<random::gpu::rand48>;

} // namespace halmd
