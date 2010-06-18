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
    template <unsigned dimension, typename sequence_type>
    void normal(sequence_type& g_val, float variance);

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
 * Fill sequence with normally distributed random numbers of given variance,
 * the value_type of sequence may be a vector
 */

template <unsigned dimension, typename sequence_type>
void random::normal(sequence_type& g_val, float variance)
{
    typedef typename sequence_type::value_type value_type;

    try {
        const uint stride = sizeof(value_type) / sizeof(float);

        // successively assign the dimensional components,
        // the dimension can not be derived from value_type
        // being often an aligned type like float2 or float4
        for (uint i=0; i < dimension; i++) {
            rng_.normal(reinterpret_cast<float*>((value_type*)g_val) + i,
                        g_val.size(), variance, stride);
        }
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to generate normal distribution on GPU");
    }
}

}} // namespace rng::gpu

} // namespace halmd

#endif /* ! HALMD_RNG_GPU_RANDOM_HPP */
