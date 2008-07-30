/* Parallel radix sort
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
#include <deque>
#include <iomanip>
#include <iostream>
#include <libgen.h>
#include <stdexcept>
#include <stdio.h>
#include <vector>

#include <radix.hpp>
#include <rand48.hpp>
#include <timer.hpp>

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
	    ("count,N", po::value<uint>(&count)->default_value(10000),
	     "number of elements")
	    ("device,D", po::value<unsigned short>(&device)->default_value(0),
	     "CUDA device")
	    ("blocks,B", po::value<uint>(&blocks)->default_value(16),
	     "number of blocks in grid")
	    ("threads,T", po::value<uint>(&threads)->default_value(128),
	     "number of threads per block")
	    ("seed,S", po::value<uint>(&seed)->default_value(42),
	     "random number generator seed")
	    ("verbose,v", po::bool_switch(&verbose),
	     "print results")
	    ("help,h", "display this help and exit");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);

	if (count < 2) {
	    throw std::logic_error("number of elements must be greater than 1");
	}
	if (threads == 0 || threads % 16) {
	    throw std::logic_error("number of threads must be a multiple of half-warp");
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
	cuda::event start, stop;

	// generate array of random integers in [0, 2^32-1] on GPU
	cuda::vector<uint> g_array(count);
	cuda::host::vector<uint> h_array(count);
	cuda::config dim((count + threads - 1) / threads, threads);
	mdsim::rand48 rng(dim);
	rng.set(seed, stream);
	rng.get(g_array, stream);
	cuda::copy(g_array, h_array, stream);
	stream.synchronize();

	// parallel radix sort
	mdsim::radix_sort<uint> radix(count, blocks, threads);
	cuda::vector<uint> g_dummy(count);
	cuda::host::vector<uint> h_array2(count);
	start.record(stream);
	radix(g_array, g_dummy, stream);
	stop.record(stream);
	cuda::copy(g_array, h_array2, stream);

	// serial radix sort
	std::vector<uint> h_array3(count);
	boost::array<std::deque<uint>, radix::BUCKET_SIZE> h_buckets;
	std::copy(h_array.begin(), h_array.end(), h_array3.begin());
	mdsim::real_timer timer;
	timer.start();
	for (uint r = 0; r < 32; r += radix::RADIX) {
	    for (uint i = 0; i < count; ++i) {
		h_buckets[(h_array3[i] >> r) & 0xff].push_back(h_array3[i]);
	    }
	    h_array3.clear();
	    for (uint b = 0; b < h_buckets.size(); ++b) {
		while (h_buckets[b].size()) {
		    h_array3.push_back(h_buckets[b].front());
		    h_buckets[b].pop_front();
		}
	    }
	}
	timer.stop();

	// wait for GPU to finish
	stream.synchronize();

	if (verbose) {
	    // write results to stdout
	    for (uint i = 0; i < count; ++i) {
		std::cout << "a[" << std::setw(6) << i << "] = " << std::setw(6)
		          << h_array[i] << ",\t[GPU] " << std::setw(10) << h_array2[i]
			  << ",\t[CPU] " << std::setw(10) << h_array3[i]
			  << ((h_array2[i] != h_array3[i]) ? " << MISMATCH\n" : "\n");
	    }
	}

	std::cout << "GPU time: " << std::fixed << std::setprecision(3)
	          << (stop - start) * 1e3 << " ms\n"
		  << "CPU time: " << std::fixed << std::setprecision(3)
		  << timer.elapsed() * 1e3 << " ms\n";

	// verify results
	if (!std::equal(h_array2.begin(), h_array2.end(), h_array3.begin())) {
	    throw std::logic_error("GPU and CPU element mismatch");
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
