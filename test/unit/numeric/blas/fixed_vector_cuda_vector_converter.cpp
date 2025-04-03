/*
 * Copyright © 2025  Felix Höfling
 * Copyright © 2012  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE fixed_vector_cuda_vector_converter
#include <boost/test/unit_test.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/version.hpp>
#include <limits>
#include <stdint.h> // uint64_t

#include <halmd/config.hpp> // HALMD_GPU_DOUBLE_PRECISION
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/cast.hpp>
#include <halmd/utility/demangle.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/unit/numeric/blas/fixed_vector_cuda_vector_converter_kernel.hpp>

using namespace std;
using namespace halmd;

#if BOOST_VERSION < 106500
BOOST_GLOBAL_FIXTURE(set_cuda_device);
#else
BOOST_TEST_GLOBAL_FIXTURE(set_cuda_device);
#endif

template <typename T>
struct epsilon;

/**
 * Conversion for single-precision floating-point is exact
 */
template <>
struct epsilon<float>
{
    double operator()() const
    {
        return numeric_limits<float>::epsilon();
    }
};

/**
 * Conversion for double-precision floating-point is lossy, since 2× float
 * has 2×24 = 48 significant bits, while double has 53 significant bits
 */
template <>
struct epsilon<double>
{
    double operator()() const
    {
        return 10 * numeric_limits<double>::epsilon();
    }
};

/**
 * Double-single precision has 2×24 = 48 significant bits
 */
template <>
struct epsilon<dsfloat>
{
    double operator()() const
    {
        return 10 * numeric_limits<double>::epsilon();
    }
};

//
// Input values
//

template <typename T>
struct index_to_value
{
    double operator()(size_t i) const
    {
        return (i + 1) * M_PI;
    }

    double operator()(size_t i, size_t j) const
    {
        return (i + 1) * M_PI + (j + 1) * M_E;
    }
};

template <>
struct index_to_value<int>
{
    int operator()(size_t i) const
    {
        return checked_narrowing_cast<int>(1 - (int64_t(i + 1) << 25));
    }
};

template <>
struct index_to_value<unsigned int>
{
    unsigned int operator()(size_t i) const
    {
        return checked_narrowing_cast<int>((uint64_t(i + 1) << 25) - 1);
    }
};

//
// Single-precision tests
//

typedef boost::mpl::vector<
    pair<fixed_vector<float, 3>, int>,
    pair<fixed_vector<float, 2>, int>,
    pair<fixed_vector<float, 3>, unsigned int>,
    pair<fixed_vector<float, 2>, unsigned int>
>::type float_int_types;

// basic conversion of a single array entry on the host, via CPU registers
BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_float_int_basic, pair_type, float_int_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();

    BOOST_TEST_MESSAGE("Splitting and merging of pair of " << demangled_name<vector_type>() << " and "
        << demangled_name<scalar_type>() << " to float4, tolerance: " << vector_eps);

    vector_type u;
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        u[j] = index_to_value<vector_value_type>()(0, j);
    }
    scalar_type v = index_to_value<scalar_type>()(0);
    float4 input;
    input <<= halmd::tie(u, v);

    volatile float w = input.w;          // make sure that 'input' is not optimized out,
    float4 output = input;               // note that tie() is not defined for 'volatile float4'
    output.w = w;

    halmd::tie(u, v) <<= output;
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(0, j), vector_eps );
    }
    BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(0) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_float_int_converter_one, pair_type, float_int_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();

    // set up input array in pinned host memory
    cuda::config dim(2, 32);
    cuda::memory::host::vector<float4> h_input(dim.threads());
    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            u[j] = index_to_value<vector_value_type>()(i, j);
        }
        scalar_type v = index_to_value<scalar_type>()(i);
        h_input[i] <<= tie(u, v);
    }

    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= h_input[i];
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(i) );
    }

    // copy to GPU memory
    cuda::memory::device::vector<float4> g_input(h_input.size());
    cuda::copy(h_input.begin(), h_input.end(), g_input.begin());
    cuda::memory::device::vector<float4> g_output(h_input.size());
    cuda::memset(g_output.begin(), g_output.end(), 0);
    cuda::memory::device::vector<vector_type> g_u(h_input.size());
    cuda::memset(g_u.begin(), g_u.end(), 0);
    cuda::memory::device::vector<scalar_type> g_v(h_input.size());
    cuda::memset(g_v.begin(), g_v.end(), 0);

    // call CUDA kernel for the conversion
    float_kernel<vector_type, scalar_type>::kernel.converter_one.configure(
        dim.grid, dim.block);
    float_kernel<vector_type, scalar_type>::kernel.converter_one(g_input,
        g_output, g_u, g_v);
    cuda::thread::synchronize();

    // copy output back to host memory
    cuda::memory::host::vector<float4> h_output(h_input.size());
    cuda::copy(g_output.begin(), g_output.end(), h_output.begin());
    cuda::memory::host::vector<vector_type> h_u(h_input.size());
    cuda::copy(g_u.begin(), g_u.end(), h_u.begin());
    cuda::memory::host::vector<scalar_type> h_v(h_input.size());
    cuda::copy(g_v.begin(), g_v.end(), h_v.begin());

    // verify result
    for (size_t i = 0; i < h_output.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= h_output[i];
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(i) );
    }

    for (size_t i = 0; i < h_u.size(); ++i) {
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( h_u[i][j], index_to_value<vector_value_type>()(i, j), vector_eps );
        }
    }

    for (size_t i = 0; i < h_v.size(); ++i) {
        BOOST_CHECK_EQUAL( h_v[i], index_to_value<scalar_type>()(i) );
    }
}

typedef boost::mpl::vector<
    pair<fixed_vector<float, 3>, float>,
    pair<fixed_vector<float, 2>, float>
>::type float_float_types;

// basic conversion of a single array entry on the host, via CPU registers
BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_float_float_basic, pair_type, float_float_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();
    double scalar_eps = epsilon<scalar_type>()();

    BOOST_TEST_MESSAGE("Splitting and merging of pair of " << demangled_name<vector_type>() << " and "
        << demangled_name<scalar_type>() << " to float4, tolerance: " << vector_eps);

    vector_type u;
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        u[j] = index_to_value<vector_value_type>()(0, j);
    }
    scalar_type v = index_to_value<scalar_type>()(0);
    float4 input;
    input <<= halmd::tie(u, v);

    volatile float w = input.w;          // make sure that 'input' is not optimized out,
    float4 output = input;               // note that tie() is not defined for 'volatile float4'
    output.w = w;

    halmd::tie(u, v) <<= output;
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(0, j), vector_eps );
    }
    BOOST_CHECK_CLOSE_FRACTION( v, index_to_value<scalar_type>()(0), scalar_eps );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_float_float_converter_one, pair_type, float_float_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();
    double scalar_eps = epsilon<scalar_type>()();

    cuda::config dim(2, 32);
    cuda::memory::host::vector<float4> h_input(dim.threads());
    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            u[j] = index_to_value<vector_value_type>()(i, j);
        }
        scalar_type v = index_to_value<scalar_type>()(i);
        h_input[i] <<= tie(u, v);
    }

    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= h_input[i];
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_CLOSE_FRACTION( v, index_to_value<scalar_type>()(i), scalar_eps );
    }

    cuda::memory::device::vector<float4> g_input(h_input.size());
    cuda::copy(h_input.begin(), h_input.end(), g_input.begin());
    cuda::memory::device::vector<float4> g_output(h_input.size());
    cuda::memset(g_output.begin(), g_output.end(), 0);
    cuda::memory::device::vector<vector_type> g_u(h_input.size());
    cuda::memset(g_u.begin(), g_u.end(), 0);
    cuda::memory::device::vector<scalar_type> g_v(h_input.size());
    cuda::memset(g_v.begin(), g_v.end(), 0);

    float_kernel<vector_type, scalar_type>::kernel.converter_one.configure(
        dim.grid, dim.block);
    float_kernel<vector_type, scalar_type>::kernel.converter_one(g_input,
        g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::memory::host::vector<float4> h_output(h_input.size());
    cuda::copy(g_output.begin(), g_output.end(), h_output.begin());
    cuda::memory::host::vector<vector_type> h_u(h_input.size());
    cuda::copy(g_u.begin(), g_u.end(), h_u.begin());
    cuda::memory::host::vector<scalar_type> h_v(h_input.size());
    cuda::copy(g_v.begin(), g_v.end(), h_v.begin());

    for (size_t i = 0; i < h_output.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= h_output[i];
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_CLOSE_FRACTION( v, index_to_value<scalar_type>()(i), scalar_eps );
    }

    for (size_t i = 0; i < h_u.size(); ++i) {
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( h_u[i][j], index_to_value<vector_value_type>()(i, j), vector_eps );
        }
    }

    for (size_t i = 0; i < h_v.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION( h_v[i], index_to_value<scalar_type>()(i), scalar_eps );
    }
}

//
// Double-precision tests
//

typedef boost::mpl::vector<
    pair<fixed_vector<dsfloat, 3>, int>,
    pair<fixed_vector<dsfloat, 2>, int>,
    pair<fixed_vector<dsfloat, 3>, unsigned int>,
    pair<fixed_vector<dsfloat, 2>, unsigned int>
#ifdef HALMD_GPU_DOUBLE_PRECISION
    ,
    pair<fixed_vector<double, 3>, int>,
    pair<fixed_vector<double, 2>, int>,
    pair<fixed_vector<double, 3>, unsigned int>,
    pair<fixed_vector<double, 2>, unsigned int>
#endif
>::type double_int_types;

// basic conversion of a single array entry on the host, via CPU registers
BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_double_int_basic, pair_type, double_int_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();

    BOOST_TEST_MESSAGE("Splitting and merging of pair of " << demangled_name<vector_type>() << " and "
        << demangled_name<scalar_type>() << " to float4, tolerance: " << vector_eps);

    vector_type u;
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        u[j] = index_to_value<vector_value_type>()(0, j);
    }
    scalar_type v = index_to_value<scalar_type>()(0);
    float4 input_hi, input_lo;
    tie(input_hi, input_lo) <<= halmd::tie(u, v);

    volatile float w_hi = input_hi.w;          // make sure that 'input_hi' is not optimized out,
    volatile float w_lo = input_lo.w;          // note that tie() is not defined for 'volatile float4'
    float4 output_hi = input_hi;
    float4 output_lo = input_lo;
    output_hi.w = w_hi;
    output_lo.w = w_lo;

    halmd::tie(u, v) <<= tie(output_hi, output_lo);
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(0, j), vector_eps );
    }
    BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(0) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_double_int_converter_two, pair_type, double_int_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();

    cuda::config dim(2, 32);
    cuda::memory::host::vector<float4> h_input(2 * dim.threads());

    for (size_t i = 0; i < dim.threads(); ++i) {
        vector_type u;
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            u[j] = index_to_value<vector_value_type>()(i, j);
        }
        scalar_type v = index_to_value<scalar_type>()(i);
        tie(h_input[i], h_input[i + dim.threads()]) <<= tie(u, v);
    }

    for (size_t i = 0; i < dim.threads(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_input[i], h_input[i + dim.threads()]);
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(u[j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(i) );
    }

    cuda::memory::device::vector<float4> g_input(2 * dim.threads());
    cuda::copy(h_input.begin(), h_input.end(), g_input.begin());
    cuda::memory::device::vector<float4> g_output(2 * dim.threads());
    cuda::memset(g_output.begin(), g_output.end(), 0);
    cuda::memory::device::vector<vector_type> g_u(dim.threads());
    cuda::memset(g_u.begin(), g_u.end(), 0);
    cuda::memory::device::vector<scalar_type> g_v(dim.threads());
    cuda::memset(g_v.begin(), g_v.end(), 0);

    double_kernel<vector_type, scalar_type>::kernel.converter_two.configure(
        dim.grid, dim.block);
    double_kernel<vector_type, scalar_type>::kernel.converter_two(g_input,
        g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::memory::host::vector<float4> h_output(2 * dim.threads());
    cuda::copy(g_output.begin(), g_output.end(), h_output.begin());
    cuda::memory::host::vector<vector_type> h_u(dim.threads());
    cuda::copy(g_u.begin(), g_u.end(), h_u.begin());
    cuda::memory::host::vector<scalar_type> h_v(dim.threads());
    cuda::copy(g_v.begin(), g_v.end(), h_v.begin());

    for (size_t i = 0; i < dim.threads(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_output[i], h_output[i + dim.threads()]);
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(u[j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(i) );
    }

    for (size_t i = 0; i < h_u.size(); ++i) {
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(h_u[i][j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
    }

    for (size_t i = 0; i < h_v.size(); ++i) {
        BOOST_CHECK_EQUAL( h_v[i], index_to_value<scalar_type>()(i) );
    }
}

typedef boost::mpl::vector<
    pair<fixed_vector<dsfloat, 3>, float>,
    pair<fixed_vector<dsfloat, 2>, float>,
    pair<fixed_vector<dsfloat, 3>, dsfloat>,
    pair<fixed_vector<dsfloat, 2>, dsfloat>
#ifdef HALMD_GPU_DOUBLE_PRECISION
    ,
    pair<fixed_vector<double , 3>, float>,
    pair<fixed_vector<double , 2>, float>,
    pair<fixed_vector<double , 3>, double>,
    pair<fixed_vector<double , 2>, double>
#endif
>::type double_float_types;

// basic conversion of a single array entry on the host, via CPU registers
BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_double_float_basic, pair_type, double_float_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();
    double scalar_eps = epsilon<scalar_type>()();

    BOOST_TEST_MESSAGE("Splitting and merging of pair of " << demangled_name<vector_type>() << " and "
        << demangled_name<scalar_type>() << " to float4, tolerance: " << vector_eps);

    vector_type u;
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        u[j] = index_to_value<vector_value_type>()(0, j);
    }
    scalar_type v = index_to_value<scalar_type>()(0);
    float4 input_hi, input_lo;
    tie(input_hi, input_lo) <<= halmd::tie(u, v);

    volatile float w_hi = input_hi.w;          // make sure that 'input_hi' is not optimized out,
    volatile float w_lo = input_lo.w;          // note that tie() is not defined for 'volatile float4'
    float4 output_hi = input_hi;
    float4 output_lo = input_lo;
    output_hi.w = w_hi;
    output_lo.w = w_lo;

    halmd::tie(u, v) <<= tie(output_hi, output_lo);
    for (size_t j = 0; j < vector_type::static_size; ++j) {
        BOOST_CHECK_CLOSE_FRACTION( u[j], index_to_value<vector_value_type>()(0, j), vector_eps );
    }
    BOOST_CHECK_CLOSE_FRACTION( double(v), index_to_value<scalar_type>()(0), scalar_eps );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_double_float_converter_two, pair_type, double_float_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();
    double scalar_eps = epsilon<scalar_type>()();

    cuda::config dim(2, 32);
    cuda::memory::host::vector<float4> h_input(2 * dim.threads());
    for (size_t i = 0; i < dim.threads(); ++i) {
        vector_type u;
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            u[j] = index_to_value<vector_value_type>()(i, j);
        }
        scalar_type v = index_to_value<scalar_type>()(i);
        tie(h_input[i], h_input[i + dim.threads()]) <<= tie(u, v);
    }

    for (size_t i = 0; i < dim.threads(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_input[i], h_input[i + dim.threads()]);
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(u[j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_CLOSE_FRACTION( double(v), index_to_value<scalar_type>()(i), scalar_eps );
    }

    cuda::memory::device::vector<float4> g_input(2 * dim.threads());
    cuda::copy(h_input.begin(), h_input.end(), g_input.begin());
    cuda::memory::device::vector<float4> g_output(2 * dim.threads());
    cuda::memset(g_output.begin(), g_output.end(), 0);
    cuda::memory::device::vector<vector_type> g_u(dim.threads());
    cuda::memset(g_u.begin(), g_u.end(), 0);
    cuda::memory::device::vector<scalar_type> g_v(dim.threads());
    cuda::memset(g_v.begin(), g_v.end(), 0);

    double_kernel<vector_type, scalar_type>::kernel.converter_two.configure(
        dim.grid, dim.block);
    double_kernel<vector_type, scalar_type>::kernel.converter_two(g_input,
        g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::memory::host::vector<float4> h_output(2 * dim.threads());
    cuda::copy(g_output.begin(), g_output.end(), h_output.begin());
    cuda::memory::host::vector<vector_type> h_u(dim.threads());
    cuda::copy(g_u.begin(), g_u.end(), h_u.begin());
    cuda::memory::host::vector<scalar_type> h_v(dim.threads());
    cuda::copy(g_v.begin(), g_v.end(), h_v.begin());

    for (size_t i = 0; i < dim.threads(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_output[i], h_output[i + dim.threads()]);
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(u[j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_CLOSE_FRACTION( double(v), index_to_value<scalar_type>()(i), scalar_eps );
    }

    for (size_t i = 0; i < h_u.size(); ++i) {
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(h_u[i][j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
    }

    for (size_t i = 0; i < h_v.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION( double(h_v[i]), index_to_value<scalar_type>()(i), scalar_eps );
    }
}
