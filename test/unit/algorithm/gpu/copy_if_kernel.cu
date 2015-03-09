/*
 * Copyright © 2015 Nicolas Höft
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

#include <stdint.h>

#include <halmd/algorithm/gpu/copy_if_kernel.cuh>
#include <test/unit/algorithm/gpu/copy_if_kernel.hpp>

template class halmd::algorithm::gpu::copy_if_wrapper<int*, int*, select_even<int> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int*, int*, select_odd<int> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int*, int*, select_all<int> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int*, int*, select_none<int> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int*, int*, select_prime<int> >;

template class halmd::algorithm::gpu::copy_if_wrapper<int8_t*, int8_t*, select_even<int8_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int8_t*, int8_t*, select_odd<int8_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int8_t*, int8_t*, select_all<int8_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int8_t*, int8_t*, select_none<int8_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int8_t*, int8_t*, select_prime<int8_t> >;

template class halmd::algorithm::gpu::copy_if_wrapper<int64_t*, int64_t*, select_even<int64_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int64_t*, int64_t*, select_odd<int64_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int64_t*, int64_t*, select_all<int64_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int64_t*, int64_t*, select_none<int64_t> >;
template class halmd::algorithm::gpu::copy_if_wrapper<int64_t*, int64_t*, select_prime<int64_t> >;
