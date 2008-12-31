/* Sample potential smoothing function in given range
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

#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include <exception>
#include <iostream>
#include <iomanip>
#include <libgen.h>
#include <gpu/ljfluid_base_glue.hpp>
namespace po = boost::program_options;
#define foreach BOOST_FOREACH

#define PROGRAM_NAME basename(argv[0])

int main(int argc, char **argv)
{
    // program options
    float r_cutoff, r_smooth;
    float3 range;
    unsigned int threads;
    unsigned short device;

    try {
	// parse command line options
	po::options_description opts("Program options");
	opts.add_options()
	    ("from", po::value<float>(&range.x)->default_value(1.1), "first endpoint of interval")
	    ("to", po::value<float>(&range.y)->default_value(1.14), "second endpoint of interval")
	    ("step", po::value<float>(&range.z)->default_value(0.001), "upper boundary for inteval step")
	    ("cutoff-distance", po::value<float>(&r_cutoff)->default_value(std::pow(2, 1 / 6.f)), "potential cutoff distance")
	    ("smooth-distance", po::value<float>(&r_smooth)->default_value(0.005), "potential smoothing distance")
	    ("device,D", po::value<unsigned short>(&device)->default_value(0), "CUDA device")
	    ("threads,T", po::value<unsigned int>(&threads)->default_value(128), "number of threads per block")
	    ("help,h", "display this help and exit");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);

	// print help
	if (vm.count("help")) {
	    std::cerr << "Usage: " << PROGRAM_NAME << " [OPTION]...\n\n" << opts << "\n";
	    return EXIT_SUCCESS;
	}
    }
    catch (std::exception const& e) {
	std::cerr << PROGRAM_NAME << ": " << e.what() << "\n";
	std::cerr << "Try `" << PROGRAM_NAME << " --help' for more information.\n";
	return EXIT_FAILURE;
    }

    try {
	// set CUDA device
	cuda::device::set(device);

	// copy constants to CUDA device symbols
	cuda::copy(r_cutoff, mdsim::gpu::ljfluid::r_cut);
	cuda::copy(std::pow(r_cutoff, 2), mdsim::gpu::ljfluid::rr_cut);
	cuda::copy(std::pow(r_smooth, -2), mdsim::gpu::ljfluid::rri_smooth);

	// CUDA execution dimensions
	unsigned int count = std::max(threads, (unsigned int)((range.y - range.x) / range.z));
	cuda::config dim((count + threads - 1) / threads, threads);

	// sample potential smoothing function in given range
	cuda::vector<float3> g_h(dim.threads());
	cuda::host::vector<float3> h_h(g_h.size());
	cuda::stream stream;
	cuda::configure(dim.grid, dim.block, stream);
	float2 f = make_float2(range.x, range.y);
	mdsim::gpu::ljfluid::sample_smooth_function(g_h, f);
	cuda::copy(g_h, h_h, stream);
	stream.synchronize();

	// write results to stdout
	foreach (float3& h, h_h) {
	    std::cout << std::scientific << std::setprecision(7) << h.x << "\t" << h.y << "\t" << h.z << "\n";
	}
	std::cout << "\n\n";
    }
    catch (cuda::error const& e) {
	std::cerr << PROGRAM_NAME << ": CUDA ERROR: " << e.what() << "\n";
	return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
