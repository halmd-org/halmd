/* Hilbert space-filling curve generation
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

#include <gpu/ljfluid_glue.hpp>
#include <gpu/radix_glue.hpp>
#include <gpu/scan_glue.hpp>
#include <rand48.hpp>
#include <timer.hpp>
#include <vector2d.hpp>
#include <vector3d.hpp>

namespace po = boost::program_options;
#define foreach BOOST_FOREACH

#define PROGRAM_NAME basename(argv[0])

#ifdef DIM_3D
typedef float4 T;
enum { dimension = 3 };
#else
typedef float2 T;
enum { dimension = 2 };
#endif

int main(int argc, char **argv)
{
    // program options
    uint blocks, threads, seed, sfc_level;
    float box;
    ushort device;

    try {
	// parse command line options
	po::options_description opts("Program options");
	opts.add_options()
	    ("box-length,L", po::value<float>(&box)->default_value(50.f),
	     "periodic simulation box length")
	    ("sfc-level,R", po::value<uint>(&sfc_level)->default_value(10),
	     "Hilbert code recursion depth")
	    ("device,D", po::value<ushort>(&device)->default_value(0),
	     "CUDA device")
	    ("blocks,B", po::value<uint>(&blocks)->default_value(16),
	     "number of blocks in grid")
	    ("threads,T", po::value<uint>(&threads)->default_value(128),
	     "number of threads per block")
	    ("seed,S", po::value<uint>(&seed)->default_value(42),
	     "random number generator seed")
	    ("help,h", "display this help and exit");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);

	if (threads == 0 || threads % 16) {
	    throw std::logic_error("number of threads must be a multiple of half-warp");
	}
	if (sfc_level > 10) {
	    throw std::logic_error("Hilbert code recursion too deep");
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
	boost::array<cuda::event, 4> start, stop;

	// copy device symbols to GPU
	cuda::copy(box, ljfluid::box);
	cuda::copy(sfc_level, ljfluid::sfc_level);

	// number of grid cells
	uint count = (1UL << (dimension * sfc_level));

	boost::array<cuda::vector<T>, 3> g_r;
	boost::array<cuda::host::vector<T>, 3> h_r;
	cuda::vector<T> g_array4(count);
	boost::array<cuda::vector<uint>, 2> g_sort;

	// compute lattice points on GPU
	cuda::config dim((count + threads - 1) / threads, threads);
	g_r[0].resize(count);
	g_r[0].reserve(dim.threads());
	start[0].record(stream);
	cuda::configure(dim.grid, dim.block, stream);
	ljfluid::lattice_simple(g_r[0], (1UL << sfc_level));
	stop[0].record(stream);
	h_r[0].resize(count);
	cuda::copy(g_r[0], h_r[0]);

	// seed GPU random number generator
	mdsim::rand48 rng(dim);
	start[1].record(stream);
	rng.set(seed, stream);
	stop[1].record(stream);

	// parallel radix sort
	cuda::vector<uint> g_bucket(blocks * threads * radix::BUCKETS_PER_THREAD);
	cuda::vector<uint> g_bucket2(g_bucket.size());
	cuda::config dim_scan((g_bucket.size() + 2 * radix::BUCKET_SIZE - 1) / (2 * radix::BUCKET_SIZE), radix::BUCKET_SIZE);
	cuda::vector<uint> g_array2(count), g_array3(count);
	cuda::vector<uint> g_block_sum(dim_scan.blocks_per_grid());
	cuda::vector<uint> g_block_sum2(g_block_sum.size());

	for (uint i = 0; i < g_sort.size(); ++i) {
	    g_r[i + 1].resize(g_r[0].size());
	    cuda::copy(g_r[i], g_r[i + 1]);

	    start[i + 2].record(stream);

	    if (i == 0) {
		// generate array of random integers in [0, 2^32-1] on GPU
		g_sort[0].resize(count);
		g_sort[0].reserve(dim.threads());
		rng.get(g_sort[0], stream);
	    }
	    else if (i == 1) {
		// generate 3D Hilbert space-filling curve
		g_sort[1].resize(count);
		g_sort[1].reserve(dim.threads());
		cuda::configure(dim.grid, dim.block, stream);
		ljfluid::sfc_hilbert_encode(g_r[i + 1], g_sort[1]);
	    }

	    for (uint r = 0; r < 32; r += radix::RADIX) {
		// compute partial radix counts
		cuda::configure(blocks, threads, threads * radix::BUCKETS_PER_THREAD * sizeof(uint), stream);
		radix::histogram_keys(g_sort[i], g_bucket, count, r);

		// parallel prefix sum over radix counts
		cuda::configure(dim_scan.grid, dim_scan.block, scan::boff(2 * dim_scan.threads_per_block()) * sizeof(uint), stream);
		scan::block_prefix_sum(g_bucket, g_bucket2, g_block_sum, g_bucket.size());
		cuda::configure(1, dim_scan.block, scan::boff(2 * dim_scan.threads_per_block() * sizeof(uint)), stream);
		scan::prefix_sum(g_block_sum, g_block_sum2, g_block_sum.size());
		cuda::configure(dim_scan.grid, dim_scan.block, stream);
		scan::add_block_sums(g_bucket2, g_bucket, g_block_sum2, g_bucket.size());

		// permute array
		cuda::configure(blocks, threads, threads * radix::BUCKETS_PER_THREAD * sizeof(uint), stream);
		radix::permute(g_sort[i], g_array2, g_r[i + 1], g_array4, g_bucket, count, r);
		cuda::copy(g_array2, g_sort[i], stream);
		cuda::copy(g_array4, g_r[i + 1], stream);
	    }

	    stop[i + 2].record(stream);

	    h_r[i + 1].resize(h_r[0].size());
	    cuda::copy(g_r[i + 1], h_r[i + 1], stream);
	}

	// wait for GPU to finish
	stream.synchronize();

	boost::array<std::string, 4> title = {{ "lattice", "seed", "random", "hilbert" }};

	// print measured GPU times
	for (uint j = 0; j < start.size(); ++j) {
	    std::cerr << title[j] << " GPU time:  \t" << std::setw(10) << std::fixed << std::setprecision(3) << (stop[j] - start[j]) * 1e3 << " ms" << std::endl;
	}

	// write results to stdout
	for (uint j = 0; j < h_r.size(); ++j) {
	    std::cout << "n(particle)\tx(" << title[j] << ")\ty(" << title[j] << ")\tz(" << title[j] << ")\n";
	    for (uint i = 0; i < count; ++i) {
		std::cout << std::scientific << i << "\t" << vector<float, dimension>(h_r[j][i]) << "\n";
	    }
	    std::cout << "\n\n";
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
