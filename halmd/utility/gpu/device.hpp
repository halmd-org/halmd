/*
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_UTILITY_GPU_DEVICE_HPP
#define HALMD_UTILITY_GPU_DEVICE_HPP

#include <lua.hpp>
#include <string>

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {

/**
 * The gpu::device module selects a GPU from the pool of available
 * GPUs (or the available subset of GPUs given as an option) and
 * prints driver and runtime version information.
 *
 * Further it stores a global value for the number of threads per
 * block, which gpu modules may use as a reasonable default.
 *
 * The GPU selection may be augmented with the external preload
 * library *nvlock* to consider only unused GPUs.
 * nvlock transparently emulates the CUDA exclusive mode feature
 * available for NVIDIA Tesla cards only, by locking the NVIDIA
 * device file upon context creation on a GPU.
 */
class device
{
public:
    static void luaopen(lua_State* L);

    device();
    ~device();

    static std::string nvidia_driver_version();
    static std::string compute_version();
    static std::string cuda_driver_version();
    static std::string cuda_runtime_version();
    //! validate CUDA execution configuration
    static cuda::config const& validate(cuda::config const& dim);

    static void* allocate(std::size_t bytes);
    static void deallocate(void* ptr);
    static void deallocate_all();
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_DEVICE_HPP */
