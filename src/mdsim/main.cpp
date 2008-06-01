/* Molecular Dynamics Simulation of a Lennard-Jones fluid
 *
 * Copyright (C) 2008  Peter Colberg
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

#include <boost/algorithm/string/join.hpp>
#include <cuda_wrapper.hpp>
#include <exception>
#include <iostream>
#include "log.hpp"
#include "mdsim.hpp"
#include "options.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "version.h"


int main(int argc, char **argv)
{
    mdsim::options opts;

    // parse program options
    try {
	opts.parse(argc, argv);
    }
    catch (mdsim::options::exit_exception const& e) {
	return e.status();
    }

    mdsim::log::init(opts);

    LOG(PROGRAM_NAME " (" PROGRAM_VERSION ")");
    LOG("variant: " << PROGRAM_VARIANT);
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif

    // print command line
    std::vector<std::string> cmd(argv, argv + argc);
    LOG("command line: " << boost::algorithm::join(cmd, " "));

    try {
	// set CUDA device for host context
	cuda::device::set(opts.device().value());
	LOG("CUDA device: " << cuda::device::get());

	// query CUDA device properties
	cuda::device::properties prop(cuda::device::get());
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

	// initialize molecular dynamics simulation
#ifdef DIM_3D
	mdsim::mdsim<3, vector3d<float> > sim(opts);
#else
	mdsim::mdsim<2, vector2d<float> > sim(opts);
#endif

	if (!opts.dry_run().value()) {
	    // run MD simulation
	    sim();
	}
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	LOG_WARNING(PROGRAM_NAME " aborted");
	return EXIT_FAILURE;
    }
    catch (std::exception const& e) {
	LOG_ERROR(e.what());
	LOG_WARNING(PROGRAM_NAME " aborted");
	return EXIT_FAILURE;
    }

    LOG(PROGRAM_NAME " exit");
    return EXIT_SUCCESS;
}
