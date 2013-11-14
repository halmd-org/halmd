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
#include <boost/nondet_random.hpp> // boost::random_device
#include <lua.hpp>
#include <iterator>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/random/gpu/rand48.hpp>

namespace halmd {
namespace random {
namespace gpu {

template <typename RandomNumberGenerator>
class random
{
public:
    typedef typename RandomNumberGenerator::rng_type rng_type;
    struct defaults;

    /**
     * Initialise random number generator.
     *
     * Get default seed from non-deterministic random number generator.
     * boost::random_device reads from /dev/urandom on GNU/Linux,
     * and the default cryptographic service provider on Windows.
     */
    random(
        unsigned int seed = boost::random_device()()
      , unsigned int blocks = defaults::blocks()
      , unsigned int threads = defaults::threads()
    );

    /**
     * Seed random number generator.
     */
    void seed(unsigned int seed);

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
        return rng_.dim.blocks_per_grid();
    }

    unsigned int threads()
    {
        return rng_.dim.threads_per_block();
    }

    RandomNumberGenerator const& rng()
    {
        return rng_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** pseudo-random number generator */
    RandomNumberGenerator rng_;
};

template <typename RandomNumberGenerator>
struct random<RandomNumberGenerator>::defaults
{
    static unsigned int blocks();
    static unsigned int threads();
};

/**
 * Shuffle sequence in-place
 */
template <typename RandomNumberGenerator>
template <typename Sequence>
void random<RandomNumberGenerator>::shuffle(Sequence& g_val)
{
    cuda::vector<unsigned int> g_sort_index;
    g_sort_index.resize(g_val.size());
    get(g_sort_index);
    radix_sort(g_sort_index.begin(), g_sort_index.end(), g_val.begin());
    cuda::thread::synchronize();
}

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RANDOM_HPP */
