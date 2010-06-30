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

#include <halmd/io/logger.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/rand48_kernel.cuh>

namespace halmd
{
namespace random { namespace gpu
{

enum { BLOCKS = 32 };
enum { THREADS = 32 << DEVICE_SCALE };

/**
 * Resolve module dependencies
 */
void random::depends()
{
    modules::depends<_Self, device_type>::required();
}

random::random(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // dependency injection
  , device(modules::fetch<device_type>(factory, vm))
{
    set_seed(vm);
}

void random::seed(unsigned int value)
{
    LOG("random number generator seed: " << value);

    try {
        rng_.resize(cuda::config(BLOCKS, THREADS));
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to change random number generator dimensions");
    }

    try {
        rng_.set(value);
        cuda::thread::synchronize();
        rng_.init_symbols(rand48_wrapper::a, rand48_wrapper::c, rand48_wrapper::state);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to seed random number generator");
    }
}

}} // namespace random::gpu

template class module<random::gpu::random>;

} // namespace halmd
