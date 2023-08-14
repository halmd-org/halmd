/*
 * Copyright © 2012 Peter Colberg
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

#include <halmd/algorithm/gpu/iota_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace detail {

/**
 * Fill range with sequentially increasing values.
 */
static __global__ void
iota(unsigned int* g_output, unsigned int size, unsigned int value)
{
    if (GTID < size) {
        g_output[GTID] = value + GTID;
    }
}

iota_kernel iota_kernel::kernel = {
    &detail::iota
};

} // namespace detail
} // namespace halmd
