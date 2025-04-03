/*
 * Copyright © 2011-2025 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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
#ifdef HALMD_WITH_GPU
#  include <cuda_wrapper/cuda_wrapper.hpp>
#endif

#include <halmd/utility/demangle.hpp>
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

#ifdef HALMD_WITH_GPU

/**
 * define test case template to test various CUDA types
 */
#define TEST_CASE_CUDA_TYPE_SIZE(t)     \
template <typename T, size_t N>         \
void test_ ## t();                      \
BOOST_AUTO_TEST_CASE( t )               \
{                                       \
    test_ ## t<float2, 2>();            \
    test_ ## t<float3, 2>();            \
    test_ ## t<float4, 2>();            \
    test_ ## t<float3, 3>();            \
    test_ ## t<float4, 3>();            \
    test_ ## t<float4, 4>();            \
}                                       \
template <typename T, size_t N>         \
void test_ ## t()                       \

#endif

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
    test_ ## t<15>();                   \
    test_ ## t<17>();                   \
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
    test_ ## t<float, 15>();            \
    test_ ## t<float, 17>();            \
    test_ ## t<double, 1>();            \
    test_ ## t<double, 2>();            \
    test_ ## t<double, 3>();            \
    test_ ## t<double, 4>();            \
    test_ ## t<double, 5>();            \
    test_ ## t<double, 8>();            \
    test_ ## t<double, 15>();           \
    test_ ## t<double, 17>();           \
    test_ ## t<signed char, 1>();       \
    test_ ## t<signed char, 2>();       \
    test_ ## t<signed char, 3>();       \
    test_ ## t<signed char, 4>();       \
    test_ ## t<signed char, 5>();       \
    /* larger sizes for signed char lead to overflow troubles in the tests */ \
    test_ ## t<uint64_t, 1>();          \
    test_ ## t<uint64_t, 2>();          \
    test_ ## t<uint64_t, 3>();          \
    test_ ## t<uint64_t, 4>();          \
    test_ ## t<uint64_t, 5>();          \
    test_ ## t<uint64_t, 8>();          \
    test_ ## t<uint64_t, 15>();         \
    test_ ## t<uint64_t, 17>();         \
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
    test_ ## t<float, double, 1>();             \
    test_ ## t<float, double, 2>();             \
    test_ ## t<float, double, 3>();             \
    test_ ## t<float, double, 4>();             \
    test_ ## t<double, double, 1>();            \
    test_ ## t<double, double, 2>();            \
    test_ ## t<double, double, 3>();            \
    test_ ## t<double, double, 4>();            \
    test_ ## t<unsigned int, double, 1>();      \
    test_ ## t<unsigned int, double, 2>();      \
    test_ ## t<unsigned int, double, 3>();      \
    test_ ## t<unsigned int, double, 4>();      \
    test_ ## t<unsigned int, int, 1>();         \
    test_ ## t<unsigned int, int, 2>();         \
    test_ ## t<unsigned int, int, 3>();         \
    test_ ## t<unsigned int, int, 4>();         \
}                                               \
template <typename T, typename S, size_t N>     \
void test_ ## t()                               \

BOOST_AUTO_TEST_SUITE( construction )

TEST_CASE_VECTOR_TYPE_SIZE( scalar )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x(1);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], 1);
    }
}

TEST_CASE_VECTOR_TYPE( initialiser_list )
{
    constexpr size_t N = 4;
    typedef fixed_vector<T, N> vector_type;

    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

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

#ifdef HALMD_WITH_GPU

TEST_CASE_CUDA_TYPE_SIZE( dsfloat_from_cuda )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<dsfloat, N>>()) << " from " << demangled_name<T>());

    fixed_vector<halmd::dsfloat, N> x;
    fixed_vector<float, N> hi, lo;
    for (size_t i = 0; i < N; ++i) {
        x[i] = (i + 1) + 1e-9 / (i + 1);
        tie(hi[i], lo[i]) = split(x[i]);
    }

    T cuda_hi(hi), cuda_lo(lo);   // convert to float2, ...
    dsfloat_ptr<T> ptr = { &cuda_hi, &cuda_lo };

    fixed_vector<dsfloat, N> y(ptr[0]);
    BOOST_CHECK_EQUAL(x, y);
}

#endif

BOOST_AUTO_TEST_SUITE_END() // construction

BOOST_AUTO_TEST_SUITE( blas1 )

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_1 )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -T(i + 1);                   // must be evaluated as signed expression if T is signed
    }

    T result = N * (N + 1) / 2;
    if (std::is_unsigned<T>::value) {       // handle unsigned types
        result = T(-1) - result + 1;
    }
    BOOST_CHECK_EQUAL(norm_1(x), result);   // floating-point arithmetic with integers is exact
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_2 )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    T const eps = std::is_floating_point<T>::value ? numeric_limits<T>::epsilon() : 1;
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -T(i + 1);
    }

    // for the 2-norm, the unsigned issue above is absent as, e.g., (-3U)^2 = 9
    BOOST_CHECK_CLOSE_FRACTION(norm_2(x), std::sqrt(N * (N + 1) * (2 * N + 1) / 6), eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_inf )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -T(i + 1);
    }

    T result = std::is_unsigned<T>::value ? T(-1) : N;     // handle unsigned types
    BOOST_CHECK_EQUAL(norm_inf(x), result);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_index_norm_inf )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -T(i + 1);
    }
    size_t result = std::is_unsigned<T>::value ? 0 : N - 1;     // handle unsigned types
    BOOST_CHECK_EQUAL(index_norm_inf(x), result);
}

BOOST_AUTO_TEST_SUITE_END() // blas1

BOOST_AUTO_TEST_SUITE( operators )

TEST_CASE_VECTOR_SIZE( vector_split_dsfloat )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<dsfloat, N>>()));

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
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()) << " and " << (demangled_name<fixed_vector<S, N>>()));

    fixed_vector<T, N> x;
    fixed_vector<S, N> y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1, 2) + pow(i + 0, 2);
        y[i] = 2 * (i + 1) * i;
    }
    x += y;

    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], pow(2 * i + 1, 2));
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_subtract )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()) << " and " << (demangled_name<fixed_vector<S, N>>()));

    fixed_vector<T, N> x;
    fixed_vector<S, N> y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1, 2) + pow(i + 0, 2);
        y[i] = 2 * (i + 1) * i;
    }
    x -= y;

    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], 1);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_multiply )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()) << " and " << demangled_name<S>());

    double const eps = std::is_floating_point<T>::value ? numeric_limits<T>::epsilon() : 0.5;
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = S(M_PI);
    x *= y;

    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(x[i], M_PI * (i + 1), eps);
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_assign_divide )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()) << " and " << demangled_name<S>());

    double const eps = std::is_floating_point<T>::value ? numeric_limits<T>::epsilon() : 0.5;
    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = S(M_PI);
    x /= y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_SMALL(x[i] - M_1_PI * (i + 1), 2 * eps);      // x[i] can become zero
    }
}

TEST_CASE_VECTOR_SIZE( vector_assign_modulus )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<unsigned int, N>>()));

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
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1, 2) + pow(i + 0, 2);
        y[i] = 2 * (i + 1) * i;
    }
    z = x + y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], pow(i + 1, 2) + pow(i + 0, 2));
        BOOST_CHECK_EQUAL(y[i], 2 * (i + 1) * i);
        BOOST_CHECK_EQUAL(z[i], pow(2 * i + 1, 2));
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_subtract )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = pow(i + 1, 2) + pow(i + 0, 2);
        y[i] = 2 * (i + 1) * i;
    }
    z = x - y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], pow(i + 1, 2) + pow(i + 0, 2));
        BOOST_CHECK_EQUAL(y[i], 2 * (i + 1) * i);
        BOOST_CHECK_EQUAL(z[i], 1);
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_sign )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y;
    for (size_t i = 0; i < N; ++i) {
        x[i] = 3 * (i + 1);
    }
    y = -x;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], 3 * (i + 1));
        BOOST_CHECK_EQUAL(y[i], -T(3) * (i + 1));     // bypass unsigned issue
    }
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_multiply_left )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = 3;
    z = y * x;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], i + 1);
        BOOST_CHECK_EQUAL(z[i], 3 * (i + 1));
    }
    BOOST_CHECK_EQUAL(y, 3);
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_multiply_right )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = 3;
    z = x * y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], i + 1);
        BOOST_CHECK_EQUAL(z[i], 3 * (i + 1));
    }
    BOOST_CHECK_EQUAL(y, 3);
}

TEST_CASE_VECTOR_TYPE_TYPE_SIZE( vector_divide )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = i + 1;
    }
    S y = 4;            // floating-point division by powers of 2 is exact
    z = x / y;
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], i + 1);
        BOOST_CHECK_EQUAL(z[i], T(i + 1) / 4);
    }
}

TEST_CASE_VECTOR_SIZE( vector_modulus )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<unsigned int, N>>()));

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
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()) << " and " << (demangled_name<fixed_vector<S, N>>()));

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
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x;
    for (size_t i = 0; i < N; ++i) {
        x[i] = -T(i + 1);
    }
    // for the inner product, there is no unsigned issue as, e.g., (-3U)^2 = 9
    BOOST_CHECK_EQUAL(inner_prod(x, x), N * (N + 1) * (2 * N + 1) / 6);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_prod )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y, z;
    T a = T(1) / 4;               // floating-point division by powers of 2 is exact
    for (size_t i = 0; i < N; ++i) {
        x[i] = 4 * (i + 1);
        y[i] = a * (i + 1);
    }
    z = element_prod(x, y);

    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], 4 * (i + 1));
        BOOST_CHECK_EQUAL(y[i], a * (i + 1));
        if (std::is_floating_point<T>::value) {
            BOOST_CHECK_EQUAL(z[i], pow(i + 1, 2));
        }
        else {
            BOOST_CHECK_EQUAL(z[i], 0);
        }
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_div )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = T(M_PI) * (i + 1);
        y[i] = T(M_PI) * (i + 1);
    }
    z = element_div(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(x[i], y[i]);
        BOOST_CHECK_EQUAL(z[i], 1);
    }
}

TEST_CASE_VECTOR_SIZE( vector_element_mod )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<unsigned int, N>>()));

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
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = (1UL << i) % N;
        y[i] = (i * i + 1) % N;
    }
    z = element_max(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(z[i], std::max(x[i], y[i]));
    }
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_element_min )
{
    BOOST_TEST_MESSAGE("Testing " << (demangled_name<fixed_vector<T, N>>()));

    fixed_vector<T, N> x, y, z;
    for (size_t i = 0; i < N; ++i) {
        x[i] = (1UL << i) % N;
        y[i] = (i * i + 1) % N;
    }
    z = element_min(x, y);
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(z[i], std::min(x[i], y[i]));
    }
}

BOOST_AUTO_TEST_SUITE_END() // operators
