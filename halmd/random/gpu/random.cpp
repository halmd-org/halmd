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

#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace random { namespace gpu
{

template <typename RandomNumberGenerator>
random<RandomNumberGenerator>::random(
    unsigned int seed
  , unsigned int blocks
  , unsigned int threads
  , unsigned int shuffle_threads
)
  // allocate random number generator state
  : rng(blocks, threads)
  , shuffle_threads_(shuffle_threads)
{
    LOG("random number generator seed: " << seed);
    try {
        rng.seed(seed);
        cuda::copy(rng.rng(), get_random_kernel<rng_type>().rng);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to seed random number generator");
        throw;
    }
}

/**
 * fill array with uniform random numbers in [0.0, 1.0)
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::uniform(cuda::vector<float>& g_v)
{
    try {
        cuda::configure(rng.dim.grid, rng.dim.block);
        get_random_kernel<rng_type>().uniform(g_v, g_v.size());
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to fill vector with uniform random numbers");
        throw;
    }
}

/**
 * fill array with random integers in [0, 2^32-1]
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::get(cuda::vector<unsigned int>& g_v)
{
    try {
        cuda::configure(rng.dim.grid, rng.dim.block);
        get_random_kernel<rng_type>().get(g_v, g_v.size());
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to fill vector with uniform integer random numbers");
        throw;
    }
}

/**
 * fill array with normal distributed random numbers in [0.0, 1.0)
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::normal(cuda::vector<float>& g_v, float mean, float sigma)
{
    try {
        cuda::configure(rng.dim.grid, rng.dim.block);
        get_random_kernel<rng_type>().normal(g_v, g_v.size(), mean, sigma);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to fill vector with normal random numbers");
        throw;
    }
}

template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(RandomNumberGenerator::name());
    module(L, "libhalmd")
    [
        namespace_("gpu")
        [
            namespace_("random")
            [
                class_<random, shared_ptr<random> >(class_name.c_str())
                    .def(constructor<
                         unsigned int
                       , unsigned int
                       , unsigned int
                       , unsigned int
                     >())
                    .property("blocks", &random::blocks)
                    .property("threads", &random::threads)
                    .scope
                    [
                        class_<defaults>("defaults")
                            .scope
                            [
                                def("seed", &defaults::seed)
                              , def("threads", &defaults::threads)
                              , def("blocks", &defaults::blocks)
                              , def("shuffle_threads", &defaults::shuffle_threads)
                            ]
                    ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_random_gpu_random(lua_State* L)
{
    random<rand48>::luaopen(L);
    return 0;
}

}} // namespace random::gpu

template class random::gpu::random<random::gpu::rand48>;

} // namespace halmd
