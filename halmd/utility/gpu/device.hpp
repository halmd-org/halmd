/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_UTILITY_GPU_DEVICE_HPP
#define HALMD_UTILITY_GPU_DEVICE_HPP

#include <lua.hpp>
#include <string>

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd
{
namespace utility { namespace gpu
{

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

    device(std::vector<int> devices, unsigned int threads);
    ~device();

    static std::string nvidia_driver_version();
    static std::string compute_version();
#if CUDA_VERSION >= 2020
    static std::string cuda_driver_version();
    static std::string cuda_runtime_version();
#endif

    cuda::device::properties const& properties() const
    {
        return prop_;
    }

    unsigned int threads() const
    {
        return threads_;
    }

private:
    /** number of CUDA threads per block */
    unsigned int threads_;
    /** device properties */
    cuda::device::properties prop_;
    /** selected CUDA device context */
    boost::shared_ptr<cuda::driver::context> context_;
};

}} // namespace utility::gpu

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_DEVICE_HPP */
