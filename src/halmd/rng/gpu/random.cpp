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
#include <halmd/rng/gpu/random.hpp>
#include <halmd/rng/gpu/rand48.cuh>
#include <halmd/rng/rand48.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/utility/module.hpp>

namespace halmd
{
namespace rng { namespace gpu
{

enum { BLOCKS = 32 };
enum { THREADS = 32 << DEVICE_SCALE };

/**
 * Resolve module dependencies
 */
void random::depends()
{
    modules::required<_Self, device_type>();
}

random::random(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , device(modules::fetch<device_type>(vm))
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

}} // namespace rng::gpu

template class module<rng::gpu::random>;

} // namespace halmd
