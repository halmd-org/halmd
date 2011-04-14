/*
 * Copyright © 2011  Peter Colberg
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

#define BOOST_TEST_MODULE variant
#include <boost/test/unit_test.hpp>

#include <algorithm> // std::max
#include <boost/mpl/sizeof.hpp>
#include <boost/type_traits/alignment_of.hpp>

#include <test/tools/cuda.hpp>
#include <test/unit/utility/gpu/variant_kernel.hpp>

using namespace std;

/**
   @file variant.cpp

   This file tests a simple implementation of a discriminated union or
   variant container. A variant container enables type-safe storage of a
   single value in a variable, where the value may have multiple types.
   An example for a sophisticated implementation is boost::variant.

   On the GPU, we make use of variants to define constants and textures
   for different spatial dimensions. As an example, the box length of the
   simulation box may be a two- or three-component vector, so we would use
   a __constant__ variant with the types float2 and float3.
 */

/**
   This test suite checks prerequisites needed for the variant container.

   The tests ensure that the metaprogramming tools work with CUDA (>=3.2),
   which includes boost::alignment_of and boost::mpl::sizeof_.
 */
BOOST_AUTO_TEST_SUITE( prerequisites )

// tests assume CUDA-supported host architecture (x86, x86_64)
BOOST_AUTO_TEST_SUITE( host )

// We could use Boost preprocessor macros (e.g. BOOST_PP_LIST_FOREACH),
// but forbear from doing so to keep the tests readable and preserve
// unique line numbers for each BOOST_CHECK_*.

BOOST_AUTO_TEST_CASE( boost_alignment_of )
{
    BOOST_CHECK_EQUAL( boost::alignment_of<char>::value, sizeof(char) );
    BOOST_CHECK_EQUAL( boost::alignment_of<short>::value, sizeof(short) );
    BOOST_CHECK_EQUAL( boost::alignment_of<int>::value, sizeof(int) );
    BOOST_CHECK_EQUAL( boost::alignment_of<long>::value, sizeof(long) );
    BOOST_CHECK_EQUAL( boost::alignment_of<float>::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::alignment_of<double>::value, sizeof(double) );

    BOOST_CHECK_EQUAL( boost::alignment_of<float1>::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::alignment_of<float2>::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<float3>::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::alignment_of<float4>::value, 16LU );

    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<char> >::value, sizeof(char) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<short> >::value, sizeof(short) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<int> >::value, sizeof(int) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<long> >::value, sizeof(long) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<float> >::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<double> >::value, sizeof(double) );

    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<float1> >::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<float2> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<float3> >::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4<float4> >::value, 16LU );

    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<char> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<short> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<int> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<long> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<float> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<double> >::value, 8LU );

    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<float1> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<float2> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<float3> >::value, 8LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align8<float4> >::value, 8LU );

    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<char> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<short> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<int> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<long> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<float> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<double> >::value, 16LU );

    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<float1> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<float2> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<float3> >::value, 16LU );
    BOOST_CHECK_EQUAL( boost::alignment_of<vector4_align16<float4> >::value, 16LU );
}

BOOST_AUTO_TEST_CASE( boost_mpl_sizeof_ )
{
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<char>::value, sizeof(char) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<short>::value, sizeof(short) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<int>::value, sizeof(int) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<long>::value, sizeof(long) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<float>::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<double>::value, sizeof(double) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<float1>::value, sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<float2>::value, 2 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<float3>::value, 3 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<float4>::value, 4 * sizeof(float) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<char> >::value, 4 * sizeof(char) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<short> >::value, 4 * sizeof(short) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<int> >::value, 4 * sizeof(int) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<long> >::value, 4 * sizeof(long) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<float> >::value, 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<double> >::value, 4 * sizeof(double) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<float1> >::value, 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<float2> >::value, 8 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<float3> >::value, 12 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4<float4> >::value, 16 * sizeof(float) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<char> >::value, max(4 * sizeof(char), 8LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<short> >::value, max(4 * sizeof(short), 8LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<int> >::value, max(4 * sizeof(int), 8LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<long> >::value, max(4 * sizeof(long), 8LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<float> >::value, max(4 * sizeof(float), 8LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<double> >::value, max(4 * sizeof(double), 8LU) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<float1> >::value, 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<float2> >::value, 8 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<float3> >::value, 12 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align8<float4> >::value, 16 * sizeof(float) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<char> >::value, max(4 * sizeof(char), 16LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<short> >::value, max(4 * sizeof(short), 16LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<int> >::value, max(4 * sizeof(int), 16LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<long> >::value, max(4 * sizeof(long), 16LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<float> >::value, max(4 * sizeof(float), 16LU) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<double> >::value, max(4 * sizeof(double), 16LU) );

    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<float1> >::value, 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<float2> >::value, 8 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<float3> >::value, 12 * sizeof(float) );
    BOOST_CHECK_EQUAL( boost::mpl::sizeof_<vector4_align16<float4> >::value, 16 * sizeof(float) );
}

BOOST_AUTO_TEST_SUITE_END() // host

// boost::mpl::sizeof_ compiles only with CUDA 3.2 or higher
#if CUDA_VERSION >= 3020

BOOST_FIXTURE_TEST_SUITE( gpu, set_cuda_device )

template <typename T>
size_t test_gpu_alignment_of()
{
    cuda::vector<size_t> g_result(1);
    cuda::host::vector<size_t> h_result(1);
    cuda::configure(1, 1);
    prerequisites_wrapper<T>::kernel.alignment_of(g_result);
    cuda::copy(g_result, h_result);
    return h_result.front();
}

BOOST_AUTO_TEST_CASE( boost_alignment_of )
{
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<char>(), sizeof(char) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<short>(), sizeof(short) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<int>(), sizeof(int) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<long>(), sizeof(long) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<float>(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<double>(), sizeof(double) );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<float1>(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<float2>(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<float3>(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<float4>(), 16LU );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<char> >(), sizeof(char) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<short> >(), sizeof(short) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<int> >(), sizeof(int) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<long> >(), sizeof(long) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<float> >(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<double> >(), sizeof(double) );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<float1> >(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<float2> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<float3> >(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4<float4> >(), 16LU );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<char> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<short> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<int> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<long> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<float> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<double> >(), 8LU );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<float1> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<float2> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<float3> >(), 8LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align8<float4> >(), 8LU );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<char> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<short> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<int> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<float> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<float> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<double> >(), 16LU );

    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<float1> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<float2> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<float3> >(), 16LU );
    BOOST_CHECK_EQUAL( test_gpu_alignment_of<vector4_align16<float4> >(), 16LU );
}

template <typename T>
size_t test_gpu_sizeof_()
{
    cuda::vector<size_t> g_result(1);
    cuda::host::vector<size_t> h_result(1);
    cuda::configure(1, 1);
    prerequisites_wrapper<T>::kernel.sizeof_(g_result);
    cuda::copy(g_result, h_result);
    return h_result.front();
}

BOOST_AUTO_TEST_CASE( boost_mpl_sizeof_ )
{
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<char>(), sizeof(char) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<short>(), sizeof(short) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<int>(), sizeof(int) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<long>(), sizeof(long) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<float>(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<double>(), sizeof(double) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<float1>(), sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<float2>(), 2 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<float3>(), 3 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<float4>(), 4 * sizeof(float) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<char> >(), 4 * sizeof(char) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<short> >(), 4 * sizeof(short) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<int> >(), 4 * sizeof(int) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<long> >(), 4 * sizeof(long) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<float> >(), 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<double> >(), 4 * sizeof(double) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<float1> >(), 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<float2> >(), 8 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<float3> >(), 12 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4<float4> >(), 16 * sizeof(float) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<char> >(), max(4 * sizeof(char), 8LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<short> >(), max(4 * sizeof(short), 8LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<int> >(), max(4 * sizeof(int), 8LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<long> >(), max(4 * sizeof(long), 8LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<float> >(), max(4 * sizeof(float), 8LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<double> >(), max(4 * sizeof(double), 8LU) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<float1> >(), 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<float2> >(), 8 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<float3> >(), 12 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align8<float4> >(), 16 * sizeof(float) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<char> >(), max(4 * sizeof(char), 16LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<short> >(), max(4 * sizeof(short), 16LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<int> >(), max(4 * sizeof(int), 16LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<long> >(), max(4 * sizeof(long), 16LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<float> >(), max(4 * sizeof(float), 16LU) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<double> >(), max(4 * sizeof(double), 16LU) );

    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<float1> >(), 4 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<float2> >(), 8 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<float3> >(), 12 * sizeof(float) );
    BOOST_CHECK_EQUAL( test_gpu_sizeof_<vector4_align16<float4> >(), 16 * sizeof(float) );
}

BOOST_AUTO_TEST_SUITE_END() // gpu

#endif /* CUDA_VERSION >= 3020 */

BOOST_AUTO_TEST_SUITE_END() // prerequisites
