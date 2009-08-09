/* Molecular Dynamics simulation of a Lennard-Jones fluid
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>
#include <fstream>
#include <iostream>
#include <ljgpu/mdlib.hpp>
#include <ljgpu/mdsim/hardsphere.hpp>
#ifdef WITH_CUDA
# include <ljgpu/mdsim/ljfluid_gpu_cell.hpp>
# include <ljgpu/mdsim/ljfluid_gpu_nbr.hpp>
# include <ljgpu/mdsim/ljfluid_gpu_square.hpp>
#endif
#include <ljgpu/mdsim/ljfluid_host.hpp>
#include <ljgpu/mdsim/mdsim.hpp>
#include <ljgpu/options.hpp>
#include <ljgpu/version.h>
#include <string>
using namespace ljgpu;

#define foreach BOOST_FOREACH

#ifdef WITH_CUDA
template <typename mdsim_backend>
static typename boost::enable_if<typename mdsim_backend::impl_type::impl_gpu, int>::type
_mdsim(options const& opt)
{
    // query NVIDIA driver version
#if !defined(__DEVICE_EMULATION__)
    std::stringbuf buf;
    std::ifstream ifs("/proc/driver/nvidia/version");
    ifs.exceptions(std::ifstream::failbit|std::ifstream::badbit);
    ifs.get(buf, '\n');
    ifs.close();
    std::string nvidia_version(buf.str());
    boost::algorithm::trim(nvidia_version);
    LOG(nvidia_version);
#endif /* ! __DEVICE_EMULATION__ */

#if !defined(__DEVICE_EMULATION__) && (CUDA_VERSION >= 2020)
    LOG("CUDA driver version: " << (cuda::driver::version() / 1000) << "." << (cuda::driver::version() / 10 % 10));
#endif
#if (CUDART_VERSION >= 2020)
    LOG("CUDA runtime version: " << (cuda::version() / 1000) << "." << (cuda::version() / 10 % 10));
#endif

#ifndef __DEVICE_EMULATION__
    // build list of chosen or available CUDA devices
    std::vector<int> devices;
    if (!opt["device"].empty()) {
	boost::multi_array<int, 1> d(opt["device"].as<boost::multi_array<int, 1> >());
	std::copy(d.begin(), d.end(), std::back_inserter(devices));
    }
    else {
	std::copy(boost::counting_iterator<int>(0),
		  boost::counting_iterator<int>(cuda::device::count()),
		  std::back_inserter(devices));
    }
    // choose first available CUDA device
    boost::shared_ptr<cuda::driver::context> ctx;
    foreach (int i, devices) {
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

    mdsim<mdsim_backend> md(opt);
#ifndef __DEVICE_EMULATION__
    LOG("GPU allocated global device memory: " << cuda::driver::mem::used() << " bytes");
    LOG("GPU available global device memory: " << cuda::driver::mem::free() << " bytes");
    LOG("GPU total global device memory: " << cuda::driver::mem::total() << " bytes");
#endif /* ! __DEVICE_EMULATION__ */

    if (opt["dry-run"].as<bool>()) {
	return LJGPU_EXIT_SUCCESS;
    }
    if (opt["daemon"].as<bool>()) {
	if (daemon(0, 0)) {
	    char const* str = strerror(errno);
	    throw std::runtime_error(std::string("failed to detach process: ") + str);
	}
    }

    return md();
}
#endif /* WITH_CUDA */

template <typename mdsim_backend>
static typename boost::disable_if<typename mdsim_backend::impl_type::impl_gpu, int>::type
_mdsim(options const& opt)
{
    // parse processor info
    std::ifstream ifs("/proc/cpuinfo");
    ifs.exceptions(std::ifstream::failbit|std::ifstream::badbit);
    std::stringbuf buf;
    ifs.get(buf, '\0');
    ifs.close();
    std::string str(buf.str());
    boost::algorithm::trim(str);
    std::vector<std::string> cpuinfo;
    boost::algorithm::split(cpuinfo, str, boost::algorithm::is_any_of("\n"));
    typedef boost::unordered_map<std::string, std::string> cpu_map;
    std::vector<cpu_map> cpus;
    foreach (std::string const& line, cpuinfo) {
	size_t pos = line.find(':');
	if (pos != std::string::npos) {
	    std::pair<std::string, std::string> tokens(line.substr(0, pos), line.substr(pos + 1));
	    boost::algorithm::trim(tokens.first);
	    boost::algorithm::trim(tokens.second);
	    if (cpus.empty() || cpus.back().find(tokens.first) != cpus.back().end()) {
		cpus.push_back(cpu_map());
	    }
	    cpus.back().insert(tokens);
	}
    }
    foreach (cpu_map& cpu, cpus) {
	LOG("CPU: " << cpu["processor"]);
	LOG("CPU family: " << cpu["cpu family"] << "  model: " << cpu["model"] << "  stepping: " << cpu["stepping"]);
	LOG("CPU model name: " << cpu["model name"]);
	LOG("CPU clock rate: " << cpu["cpu MHz"] << " MHz");
    }

    mdsim<mdsim_backend> md(opt);

    if (opt["dry-run"].as<bool>()) {
	return LJGPU_EXIT_SUCCESS;
    }
    if (opt["daemon"].as<bool>()) {
	if (daemon(0, 0)) {
	    char const* str = strerror(errno);
	    throw std::runtime_error(std::string("failed to detach process: ") + str);
	}
    }

    return md();
}

extern "C" int mdlib_mdsim(options const& opt)
{
    int const dimension = opt["dimension"].as<int>();
    if (dimension == 3) {
	return _mdsim<MDSIM_CLASS<MDSIM_IMPL, 3> >(opt);
    }
    else if (dimension == 2) {
	return _mdsim<MDSIM_CLASS<MDSIM_IMPL, 2> >(opt);
    }
    else {
	throw std::logic_error("invalid dimension: " + boost::lexical_cast<std::string>(dimension));
    }
}

extern "C" boost::program_options::options_description mdlib_options()
{
    return options::description<MDSIM_IMPL>();
}

extern "C" std::string mdlib_backend()
{
    return MDSIM_BACKEND;
}

extern "C" std::string mdlib_variant()
{
    return PROGRAM_VARIANT;
}

extern "C" std::string mdlib_version()
{
    return PROGRAM_VERSION;
}
