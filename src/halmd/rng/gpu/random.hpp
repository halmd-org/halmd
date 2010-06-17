/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_RNG_GPU_RANDOM_HPP
#define HALMD_RNG_GPU_RANDOM_HPP

#include <algorithm>
#include <iterator>

#include <boost/shared_ptr.hpp>

#include <halmd/algorithm/radix_sort.hpp>
#include <halmd/rng/random.hpp>
#include <halmd/rng/rand48.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace rng { namespace gpu
{

class random
  : public rng::random
{
public:
    // module definitions
    typedef random _Self;
    typedef rng::random _Base;
    static void options(po::options_description& desc) {};
    static void depends();
    static void select(po::options const& vm) {}

    typedef halmd::rng::rand48 random_generator;
    typedef utility::gpu::device device_type;

    shared_ptr<device_type> device;

    random(po::options const& vm);
    virtual ~random() {}
    void seed(unsigned int value);

    template <typename sequence_type>
    void shuffle(sequence_type& g_val);
    template <typename value_type>
    void normal(value_type& r1, value_type& r2, value_type sigma2);

protected:
    /** pseudo-random number generator */
    random_generator rng_;
    /** GPU radix sort */
    radix_sort<float4> radix_sort_;
};

/**
 * Shuffle sequence in-place
 */
template <typename sequence_type>
void random::shuffle(sequence_type& g_val)
{
    cuda::vector<unsigned int> g_sort_index;
    // allocate device memory
    try {
        g_sort_index.resize(g_val.size());
        g_sort_index.reserve(device->threads());

        radix_sort_.resize(g_val.size(), device->threads());
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to allocate global device memory in random::shuffle");
    }

    try {
        rng_.get(g_sort_index);
        radix_sort_(g_sort_index, g_val);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to shuffle sequence on GPU");
    }

}

/**
 * Generate two random numbers from normal distribution
 *
 * The Box-Muller transformation for generating random numbers
 * in the normal distribution was originally described in
 *
 *   G.E.P. Box and M.E. Muller, A Note on the Generation of
 *   Random Normal Deviates, The Annals of Mathematical Statistics,
 *   1958, 29, p. 610-611
 *
 * Here, we use instead the faster polar method of the Box-Muller
 * transformation, see
 *
 *   D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
 *   Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 122
 */
// template <typename value_type>
// void random::normal(value_type& x, value_type& y, value_type sigma)
// {
//     boost::uniform_01<random_generator&> variate(rng_);
//     value_type s;
//     do {
//         x = 2. * variate() - 1.;
//         y = 2. * variate() - 1.;
//         s = x * x + y * y;
//     } while (s >= 1.);
//
//     s = sigma * std::sqrt(-2. * std::log(s) / s);
//     x *= s;
//     y *= s;
// }

}} // namespace rng::gpu

} // namespace halmd

#endif /* ! HALMD_RNG_GPU_RANDOM_HPP */
