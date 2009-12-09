/* Hilbert space-filling curve generation
 *
 * Copyright (C) 2008  Peter Colberg
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
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <deque>
#include <halmd/algorithm/radix_sort.hpp>
#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/gpu/hilbert.hpp>
#include <halmd/mdsim/gpu/lattice.hpp>
#include <halmd/rng/rand48.hpp>
#include <halmd/util/timer.hpp>
#include <iomanip>
#include <iostream>
#include <libgen.h>
#include <stdexcept>
#include <stdio.h>
#include <vector>
using namespace halmd;

namespace po = boost::program_options;

#define PROGRAM_NAME basename(argv[0])

#ifdef DIM_3D
enum { dimension = 3 };
#else
enum { dimension = 2 };
#endif

int main(int argc, char **argv)
{
    // program options
    uint threads, seed, depth;
    float box;
    ushort device;

    try {
        // parse command line options
        po::options_description opts("Program options");
        opts.add_options()
            ("box-length,L", po::value<float>(&box)->default_value(50.f),
             "periodic simulation box length")
            ("depth,R", po::value<uint>(&depth)->default_value(5),
             "Hilbert code recursion depth")
            ("device,D", po::value<ushort>(&device)->default_value(0),
             "CUDA device")
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
        if (depth > 10) {
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
        // set CUDA device
        cuda::device::set(device);
        // asynchroneous GPU operations
        cuda::stream stream;
        boost::array<cuda::event, 4> start, stop;

        // copy device symbols to GPU
        cuda::copy(box, gpu::hilbert<dimension>::box);
        cuda::copy(depth, gpu::hilbert<dimension>::depth);

        // number of grid cells
        uint count = (1UL << (dimension * depth));

        boost::array<cuda::vector<float4>, 3> g_r;
        boost::array<cuda::host::vector<float4>, 3> h_r;
        boost::array<cuda::vector<uint>, 2> g_sort;

        // compute lattice points on GPU
        cuda::config dim((count + threads - 1) / threads, threads);
        g_r[0].resize(count);
        g_r[0].reserve(dim.threads());
        start[0].record(stream);
        cuda::configure(dim.grid, dim.block, stream);
        gpu::lattice<dimension>::sc(g_r[0], (1UL << depth), box);
        stop[0].record(stream);
        h_r[0].resize(count);
        cuda::copy(g_r[0], h_r[0]);

        // seed GPU random number generator
        rand48 rng(dim);
        start[1].record(stream);
        rng.set(seed, stream);
        stop[1].record(stream);

        // parallel radix sort
        radix_sort<float4> radix;
        radix.resize(count, threads);

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
                gpu::hilbert<dimension>::curve(g_r[i + 1], g_sort[1]);
            }
            // radix sort integers and particle positions
            radix(g_sort[i], g_r[i + 1], stream);
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
