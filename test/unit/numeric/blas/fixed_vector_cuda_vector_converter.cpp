/*
 * Copyright © 2012  Peter Colberg
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

#define BOOST_TEST_MODULE fixed_vector_cuda_vector_converter
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>
#include <limits>
#include <stdint.h> // uint64_t

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/cast.hpp>
#include <test/tools/ctest.hpp>
#include <test/unit/numeric/blas/fixed_vector_cuda_vector_converter_kernel.hpp>

using namespace boost;
using namespace std;
using namespace halmd;

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

typedef boost::mpl::vector<                      pair<fixed_vector<float, 3>, int         > >       float_int_types_1;
typedef boost::mpl::push_back<float_int_types_1, pair<fixed_vector<float, 2>, int         > >::type float_int_types_2;
typedef boost::mpl::push_back<float_int_types_2, pair<fixed_vector<float, 2>, unsigned int> >::type float_int_types_3;
typedef boost::mpl::push_back<float_int_types_3, pair<fixed_vector<float, 2>, unsigned int> >::type float_int_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_float_int_converter_one, pair_type, float_int_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();

    cuda::config dim(2, 32);
    cuda::host::vector<float4> h_input(dim.threads());
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

    cuda::vector<float4> g_input(h_input.size());
    cuda::copy(h_input, g_input);
    cuda::vector<float4> g_output(h_input.size());
    cuda::memset(g_output, 0);
    cuda::vector<vector_type> g_u(h_input.size());
    cuda::memset(g_u, 0);
    cuda::vector<scalar_type> g_v(h_input.size());
    cuda::memset(g_v, 0);

    cuda::configure(dim.grid, dim.block);
    float_kernel<vector_type, scalar_type>::kernel.converter_one(g_input, g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::host::vector<float4> h_output(h_input.size());
    cuda::copy(g_output, h_output);
    cuda::host::vector<vector_type> h_u(h_input.size());
    cuda::copy(g_u, h_u);
    cuda::host::vector<scalar_type> h_v(h_input.size());
    cuda::copy(g_v, h_v);

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

typedef boost::mpl::vector<                        pair<fixed_vector<float, 3>, float> >       float_float_types_1;
typedef boost::mpl::push_back<float_float_types_1, pair<fixed_vector<float, 2>, float> >::type float_float_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_float_float_converter_one, pair_type, float_float_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();
    double scalar_eps = epsilon<scalar_type>()();

    cuda::config dim(2, 32);
    cuda::host::vector<float4> h_input(dim.threads());
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

    cuda::vector<float4> g_input(h_input.size());
    cuda::copy(h_input, g_input);
    cuda::vector<float4> g_output(h_input.size());
    cuda::memset(g_output, 0);
    cuda::vector<vector_type> g_u(h_input.size());
    cuda::memset(g_u, 0);
    cuda::vector<scalar_type> g_v(h_input.size());
    cuda::memset(g_v, 0);

    cuda::configure(dim.grid, dim.block);
    float_kernel<vector_type, scalar_type>::kernel.converter_one(g_input, g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::host::vector<float4> h_output(h_input.size());
    cuda::copy(g_output, h_output);
    cuda::host::vector<vector_type> h_u(h_input.size());
    cuda::copy(g_u, h_u);
    cuda::host::vector<scalar_type> h_v(h_input.size());
    cuda::copy(g_v, h_v);

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

typedef boost::mpl::vector<                       pair<fixed_vector<double , 3>, int         > >       double_int_types_1;
typedef boost::mpl::push_back<double_int_types_1, pair<fixed_vector<double , 2>, int         > >::type double_int_types_2;
typedef boost::mpl::push_back<double_int_types_2, pair<fixed_vector<double , 3>, unsigned int> >::type double_int_types_3;
typedef boost::mpl::push_back<double_int_types_3, pair<fixed_vector<double , 2>, unsigned int> >::type double_int_types_4;
typedef boost::mpl::push_back<double_int_types_4, pair<fixed_vector<dsfloat, 3>, int         > >::type double_int_types_5;
typedef boost::mpl::push_back<double_int_types_5, pair<fixed_vector<dsfloat, 2>, int         > >::type double_int_types_6;
typedef boost::mpl::push_back<double_int_types_6, pair<fixed_vector<dsfloat, 3>, unsigned int> >::type double_int_types_7;
typedef boost::mpl::push_back<double_int_types_7, pair<fixed_vector<dsfloat, 2>, unsigned int> >::type double_int_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_double_int_converter_two, pair_type, double_int_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();

    cuda::config dim(2, 32);
    cuda::host::vector<float4> h_input;
    h_input.reserve(2 * dim.threads());
    h_input.resize(dim.threads());
    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            u[j] = index_to_value<vector_value_type>()(i, j);
        }
        scalar_type v = index_to_value<scalar_type>()(i);
        tie(h_input[i], h_input[i + h_input.size()]) <<= tie(u, v);
    }

    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_input[i], h_input[i + h_input.size()]);
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(u[j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_EQUAL( v, index_to_value<scalar_type>()(i) );
    }

    cuda::vector<float4> g_input;
    g_input.reserve(h_input.capacity());
    g_input.resize(h_input.size());
    cuda::copy(h_input, g_input, h_input.capacity());
    cuda::vector<float4> g_output;
    g_output.reserve(h_input.capacity());
    g_output.resize(h_input.size());
    cuda::memset(g_output, 0, g_output.capacity());
    cuda::vector<vector_type> g_u(h_input.size());
    cuda::memset(g_u, 0);
    cuda::vector<scalar_type> g_v(h_input.size());
    cuda::memset(g_v, 0);

    cuda::configure(dim.grid, dim.block);
    double_kernel<vector_type, scalar_type>::kernel.converter_two(g_input, g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::host::vector<float4> h_output;
    h_output.reserve(h_input.capacity());
    h_output.resize(h_input.size());
    cuda::copy(g_output, h_output, g_output.capacity());
    cuda::host::vector<vector_type> h_u(h_input.size());
    cuda::copy(g_u, h_u);
    cuda::host::vector<scalar_type> h_v(h_input.size());
    cuda::copy(g_v, h_v);

    for (size_t i = 0; i < h_output.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_output[i], h_output[i + h_output.size()]);
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

typedef boost::mpl::vector<                         pair<fixed_vector<double , 3>, float  > >       double_float_types_1;
typedef boost::mpl::push_back<double_float_types_1, pair<fixed_vector<double , 2>, float  > >::type double_float_types_2;
typedef boost::mpl::push_back<double_float_types_2, pair<fixed_vector<double , 3>, double > >::type double_float_types_3;
typedef boost::mpl::push_back<double_float_types_3, pair<fixed_vector<double , 2>, double > >::type double_float_types_4;
typedef boost::mpl::push_back<double_float_types_4, pair<fixed_vector<dsfloat, 3>, float  > >::type double_float_types_5;
typedef boost::mpl::push_back<double_float_types_5, pair<fixed_vector<dsfloat, 2>, float  > >::type double_float_types_6;
typedef boost::mpl::push_back<double_float_types_6, pair<fixed_vector<dsfloat, 3>, dsfloat> >::type double_float_types_7;
typedef boost::mpl::push_back<double_float_types_7, pair<fixed_vector<dsfloat, 2>, dsfloat> >::type double_float_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( fixed_vector_double_float_converter_two, pair_type, double_float_types )
{
    typedef typename pair_type::first_type vector_type;
    typedef typename vector_type::value_type vector_value_type;
    typedef typename pair_type::second_type scalar_type;

    double vector_eps = epsilon<vector_value_type>()();
    double scalar_eps = epsilon<scalar_type>()();

    cuda::config dim(2, 32);
    cuda::host::vector<float4> h_input;
    h_input.reserve(2 * dim.threads());
    h_input.resize(dim.threads());
    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            u[j] = index_to_value<vector_value_type>()(i, j);
        }
        scalar_type v = index_to_value<scalar_type>()(i);
        tie(h_input[i], h_input[i + h_input.size()]) <<= tie(u, v);
    }

    for (size_t i = 0; i < h_input.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_input[i], h_input[i + h_input.size()]);
        for (size_t j = 0; j < vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION( double(u[j]), index_to_value<vector_value_type>()(i, j), vector_eps );
        }
        BOOST_CHECK_CLOSE_FRACTION( double(v), index_to_value<scalar_type>()(i), scalar_eps );
    }

    cuda::vector<float4> g_input;
    g_input.reserve(h_input.capacity());
    g_input.resize(h_input.size());
    cuda::copy(h_input, g_input, h_input.capacity());
    cuda::vector<float4> g_output;
    g_output.reserve(h_input.capacity());
    g_output.resize(h_input.size());
    cuda::memset(g_output, 0, g_output.capacity());
    cuda::vector<vector_type> g_u(h_input.size());
    cuda::memset(g_u, 0);
    cuda::vector<scalar_type> g_v(h_input.size());
    cuda::memset(g_v, 0);

    cuda::configure(dim.grid, dim.block);
    double_kernel<vector_type, scalar_type>::kernel.converter_two(g_input, g_output, g_u, g_v);
    cuda::thread::synchronize();

    cuda::host::vector<float4> h_output;
    h_output.reserve(h_input.capacity());
    h_output.resize(h_input.size());
    cuda::copy(g_output, h_output, g_output.capacity());
    cuda::host::vector<vector_type> h_u(h_input.size());
    cuda::copy(g_u, h_u);
    cuda::host::vector<scalar_type> h_v(h_input.size());
    cuda::copy(g_v, h_v);

    for (size_t i = 0; i < h_output.size(); ++i) {
        vector_type u;
        scalar_type v;
        tie(u, v) <<= tie(h_output[i], h_output[i + h_output.size()]);
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
