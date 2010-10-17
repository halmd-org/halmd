/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace random { namespace gpu
{

/**
 * Assemble module options
 */
template <typename RandomNumberGenerator>
void random<RandomNumberGenerator>::options(po::options_description& desc)
{
    desc.add_options()
        ("random-blocks", po::value<unsigned int>()->default_value(default_blocks()),
         "number of CUDA blocks")
        ("random-threads", po::value<unsigned int>()->default_value(default_threads()),
         "number of CUDA threads per block")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    register_any_converter<unsigned int>();
}

template <typename RandomNumberGenerator>
random<RandomNumberGenerator>::random(
    shared_ptr<device_type> device
  , unsigned int seed
  , unsigned int blocks
  , unsigned int threads
)
  // dependency injection
  : device(device)
  // allocate random number generator state
  , rng(blocks, threads)
{
    LOG("random number generator seed: " << seed);
    try {
        rng.seed(seed);
        cuda::copy(rng.rng(), get_random_kernel<rng_type>().rng);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw runtime_error("failed to seed random number generator");
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
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw runtime_error("failed to fill vector with uniform random numbers");
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
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw runtime_error("failed to fill vector with uniform integer random numbers");
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
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw runtime_error("failed to fill vector with normal random numbers");
    }
}

template <typename Module>
static void register_lua(lua_State* L, char const* class_name)
{
    typedef typename Module::_Base _Base;
    typedef typename Module::device_type device_type;

    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("gpu")
            [
                namespace_("random")
                [
                    class_<Module, shared_ptr<_Base>, bases<_Base> >(class_name)
                        .def(constructor<shared_ptr<device_type>, unsigned int, unsigned int, unsigned int>())

                  , def("options", &Module::options)
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        bind(&register_lua<random<rand48> >, _1, "rand48")
    ];
}

}} // namespace random::gpu

template class random::gpu::random<random::gpu::rand48>;

} // namespace halmd
