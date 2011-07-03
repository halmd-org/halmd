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

#ifndef HALMD_RANDOM_GPU_RANDOM_HPP
#define HALMD_RANDOM_GPU_RANDOM_HPP

#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <iterator>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/random/gpu/rand48.hpp>

namespace halmd
{
namespace random { namespace gpu
{

template <typename RandomNumberGenerator>
class random
{
public:
    typedef typename RandomNumberGenerator::rng_type rng_type;
    struct defaults;

    RandomNumberGenerator rng; //< FIXME private?

    static void luaopen(lua_State* L);

    random(
        unsigned int seed = defaults::seed()
      , unsigned int blocks = defaults::blocks()
      , unsigned int threads = defaults::threads()
      , unsigned int shuffle_threads = defaults::shuffle_threads()
    );

    //
    // The following functions are provided for convenience.
    // Use the CUDA device functions for more flexibility.
    //
    void uniform(cuda::vector<float>& g_v);
    void get(cuda::vector<unsigned int>& g_v);
    void normal(cuda::vector<float>& g_v, float mean, float sigma);

    template <typename Sequence>
    void shuffle(Sequence& g_val);

    unsigned int blocks()
    {
        return rng.dim.blocks_per_grid();
    }

    unsigned int threads()
    {
        return rng.dim.threads_per_block();
    }

private:
    unsigned int shuffle_threads_;
};

template <typename RandomNumberGenerator>
struct random<RandomNumberGenerator>::defaults
{
    static unsigned int seed();
    static unsigned int blocks();
    static unsigned int threads();
    static unsigned int shuffle_threads();
};

/**
 * Shuffle sequence in-place
 */
template <typename RandomNumberGenerator>
template <typename Sequence>
void random<RandomNumberGenerator>::shuffle(Sequence& g_val)
{
    typedef typename Sequence::value_type value_type;
    typedef algorithm::gpu::radix_sort<value_type> sort_type;

    cuda::vector<unsigned int> g_sort_index;
    // allocate device memory
    try {
        g_sort_index.resize(g_val.size());
        g_sort_index.reserve(shuffle_threads_);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate global device memory in random::shuffle");
        throw;
    }

    sort_type sort(g_val.size(), shuffle_threads_);
    try {
        get(g_sort_index);
        sort(g_sort_index, g_val);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to shuffle sequence on GPU");
        throw;
    }
}

}} // namespace random::gpu

} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RANDOM_HPP */
