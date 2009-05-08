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

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>
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

#ifdef WITH_CUDA
template <typename mdsim_backend>
static typename boost::enable_if<typename mdsim_backend::has_gpu, void>::type
_mdsim(options const& opt)
{
    // create CUDA context and associate it with this thread
    boost::shared_ptr<cuda::context> ctx;
    try {
	// CUDA context may have been created via LD_PRELOAD
	cuda::context::device();
    }
    catch (cu::error const&) {
	ctx.reset(new cuda::context(opt["device"].as<int>()));
    }
    LOG("CUDA device: " << cuda::context::device());

    cuda::device::properties prop(cuda::context::device());
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
    LOG("GPU allocated global device memory: " << cuda::mem::used() << " bytes");
    LOG("GPU available global device memory: " << cuda::mem::free() << " bytes");
    LOG("GPU total global device memory: " << cuda::mem::total() << " bytes");

    if (opt["dry-run"].as<bool>()) {
	return;
    }
    if (opt["daemon"].as<bool>()) {
	daemon(0, 0);
    }

    md();
}
#endif /* WITH_CUDA */

template <typename mdsim_backend>
static typename boost::disable_if<typename mdsim_backend::has_gpu, void>::type
_mdsim(options const& opt)
{
    mdsim<mdsim_backend> md(opt);

    if (opt["dry-run"].as<bool>()) {
	return;
    }
    if (opt["daemon"].as<bool>()) {
	daemon(0, 0);
    }

    md();
}

extern "C" void mdlib_mdsim(options const& opt)
{
    int const dimension = opt["dimension"].as<int>();
    if (dimension == 3) {
	_mdsim<MDSIM_BACKEND<MDSIM_IMPL, 3> >(opt);
    }
    else if (dimension == 2) {
	_mdsim<MDSIM_BACKEND<MDSIM_IMPL, 2> >(opt);
    }
    else {
	throw std::logic_error("invalid dimension: " + boost::lexical_cast<std::string>(dimension));
    }
}

extern "C" options::description mdlib_options()
{
    return options_description<MDSIM_IMPL>();
}

extern "C" std::string mdlib_version()
{
    return PROGRAM_VERSION;
}
