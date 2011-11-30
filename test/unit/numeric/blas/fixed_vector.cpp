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

#define BOOST_TEST_MODULE fixed_vector
#include <boost/test/unit_test.hpp>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/blas/blas1.hpp>

using namespace boost;
using namespace std;
using namespace halmd;

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

BOOST_AUTO_TEST_SUITE( blas1 )

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_1 )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(norm_1(x), N * (N + 1) / 2, eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_2 )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(norm_2(x), sqrt(N * (N + 1) * (2 * N + 1) / 6.), eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_norm_inf )
{
    T const eps = numeric_limits<T>::epsilon();
    fixed_vector<T, N> x;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_CLOSE_FRACTION(norm_inf(x), N, eps);
}

TEST_CASE_VECTOR_TYPE_SIZE( vector_index_norm_inf )
{
    fixed_vector<T, N> x;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = -(i + 1.);
    }
    BOOST_CHECK_EQUAL(index_norm_inf(x), N - 1);
}

BOOST_AUTO_TEST_SUITE_END() // blas1
