/*
 * Copyright © 2010  Peter Colberg
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

#define BOOST_TEST_MODULE test_gpu_reduce
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

using namespace boost;
using namespace halmd;
using namespace halmd::algorithm::gpu;

/**
 * Test »32 bit integer arithmetic using double-single floating point (~44 bit).
 */
BOOST_AUTO_TEST_CASE( double_single_reduce )
{
    reduce<
        sum_            // reduce_transform
      , float           // input_type
      , float           // coalesced_input_type
      , dsfloat         // output_type
      , dsfloat         // coalesced_output_type
      , double          // host_output_type
    > sum;
    cuda::host::vector<float> h_v(make_counting_iterator(1), make_counting_iterator(12345679));
    cuda::vector<float> g_v(h_v.size());
    cuda::copy(h_v, g_v);
    BOOST_CHECK_EQUAL(int64_t(sum(g_v)), 12345678LL * 12345679 / 2);
}
