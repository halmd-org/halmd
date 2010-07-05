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

#ifndef HALMD_RANDOM_GPU_RANDOM_HPP
#define HALMD_RANDOM_GPU_RANDOM_HPP

#include <algorithm>
#include <iterator>

#include <boost/shared_ptr.hpp>

#include <halmd/algorithm/gpu/radix.hpp>
#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/random.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace random { namespace gpu
{

class random
  : public halmd::random::random
{
public:
    // module definitions
    typedef random _Self;
    typedef halmd::random::random _Base;
    static void options(po::options_description& desc) {};
    static void depends();
    static void select(po::options const& vm) {}

    typedef rand48 random_generator;
    typedef utility::gpu::device device_type;

    shared_ptr<device_type> device;

    random(modules::factory& factory, po::options const& vm);
    virtual ~random() {}
    void seed(unsigned int value);

    template <typename sequence_type>
    void shuffle(sequence_type& g_val);
    template <typename value_type>
    void normal(value_type& r1, value_type& r2, value_type sigma2);

protected:
    /** pseudo-random number generator */
    random_generator rng_;
};

/**
 * Shuffle sequence in-place
 */
template <typename sequence_type>
void random::shuffle(sequence_type& g_val)
{
    typedef typename sequence_type::value_type value_type;
    typedef algorithm::gpu::radix_sort<value_type> sort_type;

    cuda::vector<unsigned int> g_sort_index;
    // allocate device memory
    try {
        g_sort_index.resize(g_val.size());
        g_sort_index.reserve(device->threads());
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to allocate global device memory in random::shuffle");
    }

    sort_type sort(g_val.size(), device->threads());
    try {
        rng_.get(g_sort_index);
        sort(g_sort_index, g_val);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to shuffle sequence on GPU");
    }

}

}} // namespace random::gpu

} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RANDOM_HPP */
