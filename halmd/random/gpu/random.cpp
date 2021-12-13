/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2010-2011 Peter Colberg
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

#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace random {
namespace gpu {

std::shared_ptr<logger> const logger_ = std::make_shared<logger>("random (GPU)");

template <typename RandomNumberGenerator>
random<RandomNumberGenerator>::random(
    unsigned int seed
  , unsigned int blocks
  , unsigned int threads
)
  // allocate random number generator state
  : rng_(blocks, threads)
{
    LOG("random number generator type: " << rng_.name());
    LOG_DEBUG("number of CUDA execution blocks: " << blocks);
    LOG_DEBUG("number of CUDA execution threads per block: " << threads);
    random::seed(seed);
}

template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::seed(unsigned int seed)
{
    LOG("set RNG seed: " << seed);
    rng_.seed(seed);
}

/**
 * fill array with uniform random numbers in [0.0, 1.0)
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::uniform(cuda::memory::device::vector<float>& g_v)
{
    get_random_kernel<rng_type>().uniform.configure(rng_.dim.grid,
        rng_.dim.block);
    get_random_kernel<rng_type>().uniform(g_v, g_v.size(), rng_.rng());
    cuda::thread::synchronize();
}

/**
 * fill array with random integers in [0, 2^32-1]
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::get(cuda::memory::device::vector<unsigned int>& g_v)
{
    get_random_kernel<rng_type>().get.configure(rng_.dim.grid, rng_.dim.block);
    get_random_kernel<rng_type>().get(g_v, g_v.size(), rng_.rng());
    cuda::thread::synchronize();
}

/**
 * fill array with normal distributed random numbers in [0.0, 1.0)
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::normal(cuda::memory::device::vector<float>& g_v, float mean, float sigma)
{
    get_random_kernel<rng_type>().normal.configure(rng_.dim.grid,
        rng_.dim.block);
    get_random_kernel<rng_type>().normal(g_v, g_v.size(), mean, sigma,
        rng_.rng());
    cuda::thread::synchronize();
}

template <typename RandomNumberGenerator>
unsigned int random<RandomNumberGenerator>::defaults::blocks() {
    return 32;
}
template <typename RandomNumberGenerator>
unsigned int random<RandomNumberGenerator>::defaults::threads() {
    return 256;
}

template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string const class_name = RandomNumberGenerator::name();
    module(L, "libhalmd")
    [
        namespace_("random")
        [
            namespace_("gpu")
            [
                class_<random, std::shared_ptr<random> >(class_name.c_str())
                    .def(constructor<>())
                    .def(constructor<unsigned int>())
                    .def("seed", &random::seed)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_random_gpu_random(lua_State* L)
{
    random<rand48>::luaopen(L);
    return 0;
}

} // namespace gpu
} // namespace random

template class random::gpu::random<random::gpu::rand48>;

} // namespace halmd
