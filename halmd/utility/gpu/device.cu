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

#include <cub/util_allocator.cuh>
#include <cuda_wrapper/device.hpp>
#include <cuda_wrapper/error.hpp>

#include <halmd/utility/gpu/device.hpp>

namespace halmd {
namespace detail {

// instantiation of the cub allocator, the boolean parameter
// determines whether this is a object in global scope or not (it is)
// and then will skip a call to CachingDeviceAllocator::FreeAllCached()
// in the d'tor
cub::CachingDeviceAllocator caching_allocator_(true);

} // namespace detail

void* device::allocate(size_t bytes)
{
    void* ptr;
    CUDA_CALL(detail::caching_allocator_.DeviceAllocate(&ptr, bytes, cuda::device::get()));
    return ptr;
}

void device::deallocate(void* ptr)
{
    CUDA_CALL(detail::caching_allocator_.DeviceFree(ptr, cuda::device::get()));
}

void device::deallocate_all()
{
    CUDA_CALL(detail::caching_allocator_.FreeAllCached());
}

} // namespace halmd
