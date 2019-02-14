/*
 * Copyright © 2011-2012  Peter Colberg and Felix Höfling
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

#define BOOST_TEST_MODULE fixed_vector
#include <boost/test/unit_test.hpp>

#include <vector>

#include <halmd/numeric/blas/blas.hpp>
#include <test/tools/ctest.hpp>

using namespace std;
using namespace halmd;

/**
 * define test case template to test with various types
 */
#define TEST_CASE_VECTOR_TYPE(t)        \
template <typename T>                   \
void test_ ## t();                      \
BOOST_AUTO_TEST_CASE( t )               \
{                                       \
    test_ ## t<float>();                \
    test_ ## t<double>();               \
    test_ ## t<signed char>();          \
    test_ ## t<uint64_t>();             \
}                                       \
template <typename T>                   \
void test_ ## t()                       \

/**
 * define test case template to test with various sizes
 */
#define TEST_CASE_VECTOR_SIZE(t)        \
template <size_t N>                     \
void test_ ## t();                      \
BOOST_AUTO_TEST_CASE( t )               \
{                                       \
    test_ ## t<1>();                    \
    test_ ## t<2>();                    \
    test_ ## t<3>();                    \
    test_ ## t<4>();                    \
    test_ ## t<5>();                    \
    test_ ## t<8>();                    \
    test_ ## t<13>();                   \
    test_ ## t<21>();                   \
    test_ ## t<34>();                   \
    test_ ## t<55>();                   \
}                                       \
template <size_t N>                     \
void test_ ## t()                       \

/**
 * define test case template to test with various types and sizes
 */
#define TEST_CASE_VECTOR_TYPE_SIZE(t)   \
template <typename T, size_t N>         \
void test_ ## t();                      \
BOOST_AUTO_TEST_CASE( t )               \
{                                       \
    test_ ## t<float, 1>();             \
    test_ ## t<float, 2>();             \
    test_ ## t<float, 3>();             \
    test_ ## t<float, 4>();             \
    test_ ## t<float, 5>();             \
    test_ ## t<float, 8>();             \
    test_ ## t<float, 13>();            \
    test_ ## t<float, 21>();            \
    test_ ## t<float, 34>();            \
    test_ ## t<float, 55>();            \
    test_ ## t<double, 1>();            \
    test_ ## t<double, 2>();            \
    test_ ## t<double, 3>();            \
    test_ ## t<double, 4>();            \
    test_ ## t<double, 5>();            \
    test_ ## t<double, 8>();            \
    test_ ## t<double, 13>();           \
    test_ ## t<double, 21>();           \
    test_ ## t<double, 34>();           \
    test_ ## t<double, 55>();           \
}                                       \
template <typename T, size_t N>         \
void test_ ## t()                       \


/**
 * define test case template to test with various types and sizes
 */
#define TEST_CASE_VECTOR_TYPE_TYPE_SIZE(t)      \
template <typename T, typename S, size_t N>     \
void test_ ## t();                              \
BOOST_AUTO_TEST_CASE( t )                       \
{                                               \
    test_ ## t<float, float, 1>();              \
    test_ ## t<float, float, 2>();              \
    test_ ## t<float, float, 3>();              \
    test_ ## t<float, float, 4>();              \
    test_ ## t<float, float, 5>();              \
    test_ ## t<float, float, 8>();              \
    test_ ## t<float, float, 13>();             \
    test_ ## t<float, float, 21>();             \
    test_ ## t<float, float, 34>();             \
    test_ ## t<float, float, 55>();             \
    test_ ## t<float, double, 1>();             \
    test_ ## t<float, double, 2>();             \
    test_ ## t<float, double, 3>();             \
    test_ ## t<float, double, 4>();             \
    test_ ## t<float, double, 5>();             \
    test_ ## t<float, double, 8>();             \
    test_ ## t<float, double, 13>();            \
    test_ ## t<float, double, 21>();            \
    test_ ## t<float, double, 34>();            \
    test_ ## t<float, double, 55>();            \
    test_ ## t<float, double, 1>();             \
    test_ ## t<float, double, 2>();             \
    test_ ## t<float, double, 3>();             \
    test_ ## t<float, double, 4>();             \
    test_ ## t<float, double, 5>();             \
    test_ ## t<float, double, 8>();             \
    test_ ## t<float, double, 13>();            \
    test_ ## t<float, double, 21>();            \
    test_ ## t<float, double, 34>();            \
    test_ ## t<float, double, 55>();            \
    test_ ## t<double, double, 1>();            \
    test_ ## t<double, double, 2>();            \
    test_ ## t<double, double, 3>();            \
    test_ ## t<double, double, 4>();            \
    test_ ## t<double, double, 5>();            \
    test_ ## t<double, double, 8>();            \
    test_ ## t<double, double, 13>();           \
    test_ ## t<double, double, 21>();           \
    test_ ## t<double, double, 34>();           \
    test_ ## t<double, double, 55>();           \
}                                               \
template <typename T, typename S, size_t N>     \
void test_ ## t()                               \

BOOST_AUTO_TEST_SUITE( construction )

TEST_CASE_VECTOR_TYPE_SIZE( scalar )
{
    fixed_vector<T, N> x(1);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], 1);
    }
}

TEST_CASE_VECTOR_TYPE( initialiser_list )
{
    constexpr size_t N = 4;
    typedef fixed_vector<T, N> vector_type;

    // construction by initialiser list for all class members
    vector_type x{1, 2, 3, 4};
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], T(i + 1));
    }

    // assignment by initialiser list
    x = {2, 3, 5, 7};
    T y[4] = { 2, 3, 5, 7 };
    BOOST_CHECK_EQUAL_COLLECTIONS(begin(x), end(x), y, y + 4);

    // partial initialisation
    // x = { 1 };   // compiles, but triggers an assertion if compiled in Debug mode
}

BOOST_AUTO_TEST_SUITE_END() // construction

BOOST_AUTO_TEST_SUITE( blas1 )

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_1 )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(norm_1(x), N * (N + 1) / 2, eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_2 )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(norm_2(x), sqrt(N * (N + 1) * (2 * N + 1) / 6.), eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_inf )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(norm_inf(x), N, eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_index_norm_inf )
{
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_EQUAL(index_norm_inf(x), N - 1);
}

BOOST_AUTO_TEST_SUITE_END() // blas1

BOOST_AUTO_TEST_SUITE( operators )

TEST_CASE_VECTOR_SIZE( vector_split_dsfloat )
{
    float const eps = numeric_limits<float>::epsilon();
    fixed_vector<dsfloat, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = (i + 1) + 1e-9 * (i + 1);
    }
    fixed_vector<float, N> hi, lo;
    tie(hi, lo) = split(x);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(hi[i], (i + 1.f), eps);
        BOOST_CHECK_CLOSE_FRACTION(lo[i], (i + 1.f) * 1e-9, eps);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_add )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    fixed_vector<S, N> y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1., 2) + pow(i + 0., 2);
        y[i] = 2 * (i + 1) * i;
    }
    x += y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], pow(2 * i + 1., 2), eps);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_subtract )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    fixed_vector<S, N> y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1., 2) + pow(i + 0., 2);
        y[i] = 2 * (i + 1) * i;
    }
    x -= y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], 1, eps);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_multiply )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = M_PI;
    x *= y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], M_PI * (i + 1), eps);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_divide )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = M_PI;
    x /= y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], M_1_PI * (i + 1), eps);
    }
}

TEST_CASE_VECTOR_SIZE( vector_assign_modulus )
{
    fixed_vector<unsigned int, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    int y = 3;
    x %= y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], (i + 1) % 3);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_add )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1., 2) + pow(i + 0., 2);
        y[i] = 2 * (i + 1) * i;
    }
    z = x + y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], pow(i + 1., 2) + pow(i + 0., 2), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], 2 * (i + 1) * i, eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], pow(2 * i + 1., 2), eps);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_subtract )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1., 2) + pow(i + 0., 2);
        y[i] = 2 * (i + 1) * i;
    }
    z = x - y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], pow(i + 1., 2) + pow(i + 0., 2), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], 2 * (i + 1) * i, eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], 1, eps);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_sign )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = M_PI * (i + 1);
    }
    y = -x;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], M_PI * (i + 1), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], -M_PI * (i + 1), eps);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_multiply_left )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = M_PI;
    z = y * x;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], i + 1, eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], M_PI * (i + 1), eps);
    }
    BOOST_CHECK_CLOSE_FRACTION(y, M_PI, eps);
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_multiply_right )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = M_PI;
    z = x * y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], i + 1, eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], M_PI * (i + 1), eps);
    }
    BOOST_CHECK_CLOSE_FRACTION(y, M_PI, eps);
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_divide )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = M_PI;
    z = x / y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], i + 1, eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], M_1_PI * (i + 1), eps);
    }
    BOOST_CHECK_CLOSE_FRACTION(y, M_PI, eps);
}

TEST_CASE_VECTOR_SIZE( vector_modulus )
{
    fixed_vector<unsigned int, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    int y = 3;
    z = x % y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], i + 1);
        BOOST_CHECK_EQUAL(z[i], (i + 1) % 3);
    }
    BOOST_CHECK_EQUAL(y, 3);
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_comparison )
{
    T const eps = numeric_limits<T>::epsilon(); // assume that T has lower precision than S
    fixed_vector<T, N> x;
    fixed_vector<S, N> y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = 1 + i * eps;
        y[i] = 1 + i * eps;
    }
    BOOST_CHECK_EQUAL(x, y);

    x[N / 2] = 0;
    BOOST_CHECK(!(x == y));
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_inner_prod )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(inner_prod(x, x), N * (N + 1) * (2 * N + 1) / 6., eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_prod )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = M_PI * (i + 1);
        y[i] = M_1_PI * (i + 1);
    }
    z = element_prod(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], M_PI * (i + 1), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], M_1_PI * (i + 1), eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], pow(i + 1., 2), eps);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_div )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = M_PI * (i + 1);
        y[i] = M_PI * (i + 1);
    }
    z = element_div(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], M_PI * (i + 1), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], M_PI * (i + 1), eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], 1, eps);
    }
}

TEST_CASE_VECTOR_SIZE( vector_element_mod )
{
    fixed_vector<unsigned int, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 23;
        y[i] = i + 17;
    }
    z = element_mod(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], i + 23);
        BOOST_CHECK_EQUAL(y[i], i + 17);
        BOOST_CHECK_EQUAL(z[i], 6u);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_max )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = sin((i + 1) / M_PI);
        y[i] = sin((i + 1) / M_PI + M_PI);
    }
    z = element_max(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], sin((i + 1) / M_PI), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], -sin((i + 1) / M_PI), 20 * eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], abs(sin((i + 1) / M_PI)), 20 * eps);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_min )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = sin((i + 1) / M_PI);
        y[i] = sin((i + 1) / M_PI + M_PI);
    }
    z = element_min(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], sin((i + 1) / M_PI), eps);
        BOOST_CHECK_CLOSE_FRACTION(y[i], -sin((i + 1) / M_PI), 20 * eps);
        BOOST_CHECK_CLOSE_FRACTION(z[i], -abs(sin((i + 1) / M_PI)), 20 * eps);
    }
}

BOOST_AUTO_TEST_SUITE_END() // operators
