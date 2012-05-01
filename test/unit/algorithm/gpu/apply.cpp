/*
 * Copyright © 2011  Felix Höfling
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

#define BOOST_TEST_MODULE apply
#include <boost/test/unit_test.hpp>

#include <halmd/algorithm/gpu/apply_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>

using namespace halmd;
using namespace halmd::algorithm::gpu;

BOOST_GLOBAL_FIXTURE( set_cuda_device )

/**
 * Test apply wrapper for square transform.
 */
BOOST_AUTO_TEST_CASE( apply_square )
{
    typedef fixed_vector<float, 2> vector_type;

    typedef apply_wrapper<
        square_               // transform
      , vector_type           // input_type
      , float2                  // coalesced_input_type
      , float                   // output_type
      , float                   // coalesced_output_type
    > square;

    unsigned int size = 1234567;
    cuda::host::vector<float2> h_input(size);
    cuda::host::vector<float> h_output(size);
    cuda::vector<float2> g_input(size);
    cuda::vector<float> g_output(size);

    // create sequence of 2-dim vectors
    for (unsigned int i = 0; i < size; ++i) {
        h_input[i].x = 2 * i;
        h_input[i].x = 2 * i + 1;
    }
    cuda::copy(h_input, g_input);

    try {
        cuda::configure(50, 256); // 50 blocks, 256 threads
        square::kernel.apply(g_input, g_output, size);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        BOOST_ERROR("streaming of apply kernel failed");
        throw;
    }

    cuda::copy(g_output, h_output);

    for (unsigned int i = 0; i < size; ++i) {
        float result = inner_prod(static_cast<vector_type>(h_input[i]), static_cast<vector_type>(h_input[i]));
        BOOST_CHECK_MESSAGE(h_output[i] == result,
            "h_output[" << i << "] = " << h_output[i] << " does not equal " << result
        );
    }
}

