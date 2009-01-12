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

#include <H5Cpp.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/foreach.hpp>
#include <cuda_wrapper.hpp>
#include <exception>
#include <iostream>
#include <libgen.h>
#include <ljgpu/mdlib.hpp>
#include <ljgpu/options.hpp>
#include <ljgpu/util/log.hpp>
#include <ljgpu/version.h>
using namespace boost;
using namespace std;

int main(int argc, char **argv)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    // parse program options
    ljgpu::options opt;
    try {
	opt.parse(argc, argv);
    }
    catch (ljgpu::options::exit_exception const& e) {
	return e.status();
    }

    string const backend(opt["backend"].as<string>());
    ljgpu::mdlib mdlib;
    try {
	mdlib.open(dirname(argv[0]) + ("/libljgpu_" + backend) + ".so");
    }
    catch (std::exception const& e) {
	cerr << e.what() << endl;
	return EXIT_FAILURE;
    }

    // parse backend options
    try {
	opt.parse(mdlib.options());
    }
    catch (ljgpu::options::exit_exception const& e) {
	return e.status();
    }

    ljgpu::log::init(opt["output"].as<std::string>() + ".log", opt["verbose"].as<int>());

    LOG(PROGRAM_NAME " " PROGRAM_VERSION);
    LOG("variant: " << PROGRAM_VARIANT);
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif

    // print command line
    vector<string> cmd(argv, argv + argc);
    LOG("command line: " << boost::algorithm::join(cmd, " "));

    LOG("MD simulation backend: " << backend);

#ifdef NDEBUG
    try {
#endif
	// bind process to CPU core(s)
	if (!opt["processor"].empty()) {
	    cpu_set_t cpu_set;
	    CPU_ZERO(&cpu_set);
	    BOOST_FOREACH(int cpu, opt["processor"].as<vector<int> >()) {
		LOG("adding CPU core " << cpu << " to process CPU affinity mask");
		CPU_SET(cpu, &cpu_set);
	    }
	    if (0 != sched_setaffinity(getpid(), sizeof(cpu_set_t), &cpu_set)) {
		throw logic_error("failed to set process CPU affinity mask");
	    }
	}

	// set CUDA device for host context
	if (!opt["device"].empty()) {
	    int dev = opt["device"].as<int>();
	    if (opt["device"].defaulted()) {
		char const* env = getenv("CUDA_DEVICE");
		if (env != NULL && *env != '\0') {
		    char* endptr;
		    int i = strtol(env, &endptr, 10);
		    if (*endptr != '\0' || i < 0) {
			throw logic_error(string("CUDA_DEVICE environment variable invalid: ") + env);
		    }
		    dev = i;
		}
	    }
	    cuda::device::set(dev);
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
	}

	// run MD simulation
	mdlib.mdsim(opt);
#ifdef NDEBUG
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
#endif /* NDEBUG */

    LOG(PROGRAM_NAME " exit");
    return EXIT_SUCCESS;
}
