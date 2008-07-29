/* Parallel GPU rand48 random number generator test
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
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <gpu/scan_glue.hpp>
#include <rand48.hpp>
#include <iomanip>
#include <iostream>
#include <libgen.h>
#include <stdexcept>
#include <stdio.h>
#include <timer.hpp>
#include <vector>
namespace po = boost::program_options;
#define foreach BOOST_FOREACH

#define PROGRAM_NAME basename(argv[0])

int main(int argc, char **argv)
{
    // program options
    uint count, blocks, threads, seed;
    unsigned short device;
    bool verbose;

    try {
	// parse command line options
	po::options_description opts("Program options");
	opts.add_options()
	    ("count,N", po::value<uint>(&count)->default_value(10000000),
	     "random number count")
	    ("seed,S", po::value<uint>(&seed)->default_value(42),
	     "random number generator seed")
	    ("device,D", po::value<unsigned short>(&device)->default_value(0),
	     "CUDA device")
	    ("blocks,B", po::value<uint>(&blocks)->default_value(64),
	     "number of blocks in grid")
	    ("threads,T", po::value<uint>(&threads)->default_value(128),
	     "number of threads per block")
	    ("verbose,v", po::bool_switch(&verbose),
	     "print results")
	    ("help,h", "display this help and exit");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);

	if (count < 1) {
	    throw std::logic_error("number of elements must be non-zero");
	}
	if (blocks < 1) {
	    throw std::logic_error("number of blocks must be non-zero");
	}
	if (threads < 16) {
	    throw std::logic_error("number of threads must be at least 16");
	}

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
	using namespace mdsim::gpu;

	// set CUDA device
	cuda::device::set(device);
	// asynchroneous GPU operations
	cuda::stream stream;
	boost::array<cuda::event, 2> start, stop;

	std::cout << "integers:\t\t" << std::setw(13) << count << "\n";
	std::cout << "blocks: \t\t" << std::setw(13) << blocks << "\n";
	std::cout << "threads:\t\t" << std::setw(13) << threads << "\n";
	std::cout << "seed:   \t\t" << std::setw(13) << seed << "\n";
	std::cout << std::endl;

	// seed GPU random number generator
	cuda::config dim(blocks, threads);
	mdsim::rand48 rng(dim);
	start[0].record(stream);
	rng.set(seed, stream);
	stop[0].record(stream);

	// parallel GPU rand48
	cuda::vector<uint> g_array(count);
	start[1].record(stream);
	rng.get(g_array, stream);
	stop[1].record(stream);

	cuda::host::vector<uint> h_array(count);
	cuda::copy(g_array, h_array, stream);

	// serial GNU C library rand48
	std::vector<uint> h_array2(count);
	srand48(seed);
	mdsim::real_timer timer;
	timer.start();
	std::generate(h_array2.begin(), h_array2.end(), mrand48);
	timer.stop();

	// wait for GPU to finish
	stream.synchronize();

	if (verbose) {
	    // write results to stdout
	    for (uint i = 0; i < count; ++i) {
		std::cout << "[" << std::setw(6) << i << "]\t[GPU] " << std::setw(10)
		          << h_array[i] << ",\t[CPU] " << std::setw(10) << h_array2[i]
			  << ((h_array[i] != h_array2[i]) ? " << MISMATCH\n" : "\n");
	    }
	}

	// print measured GPU times
	boost::array<std::string, 4> title = {{ "seed", "rand48", }};
	for (uint j = 0; j < start.size(); ++j) {
	    std::cout << title[j] << " GPU time:  \t" << std::setw(10) << std::fixed << std::setprecision(3) << (stop[j] - start[j]) * 1e3 << " ms" << std::endl;
	}
	// print measured CPU time
	std::cout << "rand48 CPU time:  \t" << std::setw(10) << std::fixed << std::setprecision(3) << timer.elapsed() * 1e3 << " ms" << std::endl;

	// verify results
	if (!std::equal(h_array.begin(), h_array.end(), h_array2.begin())) {
	    throw std::logic_error("GPU and CPU random number mismatch");
	}
    }
    catch (cuda::error const& e) {
	std::cerr << PROGRAM_NAME << ": CUDA ERROR: " << e.what() << "\n";
	return EXIT_FAILURE;
    }
    catch (std::exception const& e) {
	std::cerr << PROGRAM_NAME << ": ERROR: " << e.what() << "\n";
	return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
