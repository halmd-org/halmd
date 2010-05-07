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

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <exception>
#include <fstream>

#include <halmd/io/logger.hpp>
#include <halmd/utility/gpu/device.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace std;

namespace halmd
{
namespace utility { namespace gpu
{

/**
 * Assemble module options
 */
void device::options(po::options_description& desc)
{
    desc.add_options()
#ifndef __DEVICE_EMULATION__
        ("device,D", po::value<boost::multi_array<int, 1> >(),
         "CUDA device(s)")
#endif
        ("threads,T", po::value<unsigned int>()->default_value(128),
         "number of CUDA threads per block")
        ;
}

/**
 * Initialize CUDA device
 */
device::device(po::options const& vm)
{
#ifndef __DEVICE_EMULATION__
    try {
        LOG("NVIDIA driver version: " << device::nvidia_driver_version());
    }
    catch (runtime_error const& e) {
        LOG_WARNING(e.what());
    }
# if CUDA_VERSION >= 2020
    LOG("CUDA driver version: " << device::cuda_driver_version());
# endif
# if (CUDART_VERSION >= 2020)
    LOG("CUDA runtime version: " << device::cuda_runtime_version());
# endif

    // build list of chosen or available CUDA devices
    vector<int> devices;
    if (!vm["device"].empty()) {
        multi_array<int, 1> v(vm["device"].as<multi_array<int, 1> >());
        copy(v.begin(), v.end(), back_inserter(devices));
    }
    else {
        copy(counting_iterator<int>(0),
             counting_iterator<int>(cuda::device::count()),
             back_inserter(devices));
    }

    // choose first available CUDA device
    shared_ptr<cuda::driver::context> ctx;
    BOOST_FOREACH (int i, devices) {
        try {
            // create CUDA context and associate it with this thread
            ctx.reset(new cuda::driver::context(i));
            break;
        }
        catch (cuda::driver::error const&) {
            // device is compute-exlusive mode and in use
        }
    }

    LOG("CUDA device: " << cuda::driver::context::device());
    cuda::device::properties prop(cuda::driver::context::device());
#else /* ! __DEVICE_EMULATION__ */
    LOG("CUDA device: " << cuda::device::get());
    cuda::device::properties prop(cuda::device::get());
#endif /* ! __DEVICE_EMULATION__ */

    LOG("CUDA device name: " << prop.name());
    LOG("CUDA device total global memory: " << prop.total_global_mem() << " bytes");
    LOG("CUDA device shared memory per block: " << prop.shared_mem_per_block() << " bytes");
    LOG("CUDA device registers per block: " << prop.regs_per_block());
    LOG("CUDA device warp size: " << prop.warp_size());
    LOG("CUDA device maximum number of threads per block: " << prop.max_threads_per_block());
    LOG("CUDA device total constant memory: " << prop.total_const_mem());
    LOG("CUDA device major revision: " << prop.major());
    LOG("CUDA device minor revision: " << prop.minor());
    LOG("CUDA device clock frequency: " << prop.clock_rate() << " kHz");

    threads_ = vm["threads"].as<unsigned int>();
    if (threads_ < 1) {
        throw runtime_error("invalid number of CUDA threads");
    }
    if (threads_ > prop.max_threads_per_block()) {
        throw runtime_error("number of CUDA threads exceeds maximum number of threads per block");
    }
    if (threads_ & (threads_ - 1)) {
        LOG_WARNING("number of CUDA threads not a power of 2");
    }
    if (threads_ % prop.warp_size()) {
        LOG_WARNING("number of CUDA threads not a multiple of warp size");
    }

    LOG("number of CUDA threads: " << threads_);
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

#if CUDA_VERSION >= 2020

/**
 * Query CUDA driver version
 */
string device::cuda_driver_version()
{
    int major = cuda::driver::version() / 1000;
    int minor = cuda::driver::version() / 10 % 10;
    return lexical_cast<string>(major) + "." + lexical_cast<string>(minor);
}

#endif /* CUDA_VERSION >= 2020 */

#if CUDART_VERSION >= 2020

/**
 * Query CUDA runtime version
 */
string device::cuda_runtime_version()
{
    int major = cuda::version() / 1000;
    int minor = cuda::version() / 10 % 10;
    return lexical_cast<string>(major) + "." + lexical_cast<string>(minor);
}

#endif /* CUDART_VERSION >= 2020 */

}} // namespace utility::gpu

template class module<utility::gpu::device>;

} // namespace halmd
