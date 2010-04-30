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

#include <string>

#include <cuda_wrapper.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd { namespace utility { namespace gpu
{

/**
 * CUDA device configuration
 */
class device
{
public:
    typedef device module_type;

    static void options(po::options_description& desc);
    static void resolve(po::options const& vm) {}
    device(po::options const& vm);
    virtual ~device() {}
    unsigned int threads() { return threads_; }

protected:
    static std::string nvidia_driver_version();
#if CUDA_VERSION >= 2020
    static std::string cuda_driver_version();
#endif
#if CUDART_VERSION >= 2020
    static std::string cuda_runtime_version();
#endif

    /** number of CUDA threads per block */
    unsigned int threads_;
};

}}} // namespace halmd::utility::gpu

#endif /* ! HALMD_UTILITY_GPU_DEVICE_HPP */
