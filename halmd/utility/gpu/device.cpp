/*
 * Copyright © 2008-2012 Peter Colberg
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

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <exception>
#include <fstream>

#include <halmd/config.hpp> // HALMD_GPU_ARCH
#include <halmd/io/logger.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace std;

namespace halmd {

/**
 * Initialize CUDA device
 */
device::device()
{
    try {
        LOG("NVIDIA driver version: " << device::nvidia_driver_version());
    }
    catch (runtime_error const& e) {
        LOG_WARNING(e.what());
    }
    LOG("CUDA driver version: " << device::cuda_driver_version());
    LOG("CUDA runtime version: " << device::cuda_runtime_version());

    // choose first available CUDA device
    cuda::thread::synchronize();

    cuda::device::properties prop(cuda::device::get());

    LOG("GPU: " << cuda::device::get());
    LOG("GPU name: " << prop.name());
    LOG("GPU total global memory: " << prop.total_global_mem() << " bytes");
    LOG("GPU shared memory per block: " << prop.shared_mem_per_block() << " bytes");
    LOG("GPU registers per block: " << prop.regs_per_block());
    LOG("GPU warp size: " << prop.warp_size());
    LOG("GPU maximum number of threads per block: " << prop.max_threads_per_block());
    LOG("GPU total constant memory: " << prop.total_const_mem());
    LOG("GPU clock frequency: " << prop.clock_rate() << " kHz");
    LOG("CUDA compute capability: " << prop.major() << "." << prop.minor());
    LOG("CUDA compute version: " << device::compute_version());
}

/**
 * Detach CUDA runtime from CUDA device context
 *
 * This explicit clean-up is needed with CUDA < 3.0.
 */
device::~device()
{
    LOG_DEBUG("detach from CUDA device context");
    cuda::thread::exit();
}

cuda::config const& device::validate(cuda::config const& dim)
{
    cuda::thread::synchronize();
    cuda::device::properties prop(cuda::device::get());
    unsigned int threads = dim.threads_per_block();

    if (threads < 1) {
        throw runtime_error("invalid number of GPU threads: " + lexical_cast<string>(threads));
    }
    if (threads > prop.max_threads_per_block()) {
        throw runtime_error("number of GPU threads exceeds maximum threads per block: " + lexical_cast<string>(threads));
    }
    if (threads % prop.warp_size()) {
        throw runtime_error("number of GPU threads not a multiple of warp size: " + lexical_cast<string>(threads));
    }
    if (threads & (threads - 1)) {
        LOG_WARNING("number of GPU threads not a power of 2: " << threads);
    }
    return dim;
}

/**
 * Query NVIDIA driver version
 */
string device::nvidia_driver_version()
{
    string s;
    try {
        stringbuf buf;
        ifstream ifs("/proc/driver/nvidia/version");
        ifs.exceptions(std::ifstream::failbit|std::ifstream::badbit);
        ifs.get(buf, '\n');
        ifs.close();
        s = buf.str();
    }
    catch (ifstream::failure&) {
        throw runtime_error("failed to query NVIDIA driver version");
    }
    size_t pos = s.find(": ");
    if (pos != string::npos) {
        s = s.substr(pos + 2);
    }
    trim(s);
    return s;
}

/**
 * Query CUDA compute version
 */
string device::compute_version()
{
    int major = HALMD_GPU_ARCH / 100;
    int minor = HALMD_GPU_ARCH / 10 % 10;
    return lexical_cast<string>(major) + "." + lexical_cast<string>(minor);
}

/**
 * Query CUDA driver version
 */
string device::cuda_driver_version()
{
    int version = cuda::driver_version();
    int major = version / 1000;
    int minor = version / 10 % 10;
    return lexical_cast<string>(major) + "." + lexical_cast<string>(minor);
}

/**
 * Query CUDA runtime version
 */
string device::cuda_runtime_version()
{
    int version = cuda::runtime_version();
    int major = version / 1000;
    int minor = version / 10 % 10;
    return lexical_cast<string>(major) + "." + lexical_cast<string>(minor);
}

/**
 * Translate CUDA exception to Lua error message
 */
static void translate_cuda_error(lua_State* L, cuda::error const& e)
{
    lua_pushliteral(L, "[CUDA] ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

void device::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            namespace_("gpu")
            [
                class_<device, boost::shared_ptr<device> >("device")
                    .def(constructor<>())
                    .scope
                    [
                        def("nvidia_driver_version", &device::nvidia_driver_version)
                      , def("compute_version", &device::compute_version)
                      , def("cuda_driver_version", &device::cuda_driver_version)
                      , def("cuda_runtime_version", &device::cuda_runtime_version)
                    ]
            ]
        ]
    ];
    register_exception_handler<cuda::error>(&translate_cuda_error);
}

HALMD_LUA_API int luaopen_libhalmd_utility_gpu_device(lua_State* L)
{
    device::luaopen(L);
    return 0;
}

} // namespace halmd
