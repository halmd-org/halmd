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

#define BOOST_TEST_MODULE multi_range
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/multi_array.hpp>
#include <halmd/algorithm/multi_range.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace std;

template <typename T, size_t N>
bool element_equal(multi_array<T, N> const& tensor, T const& value, array<size_t, N> const& index)
{
    return tensor(index) == value;
}

struct multi_array_fixture
{
    multi_array_fixture() : size(list_of(3)(7)(17)), tensor(size)
    {
        for (size_t i = 0, n = 1; i < size[0]; ++i) {
            for (size_t j = 0; j < size[1]; ++j) {
                for (size_t k = 0; k < size[2]; ++k, ++n) {
                    tensor[i][j][k] = n;
                }
            }
        }
    }

    array<size_t, 3> const size;
    multi_array<double, 3> tensor;
};

/**
 * find each element by value
 */
BOOST_FIXTURE_TEST_CASE( find_if_every_element, multi_array_fixture )
{
    array<size_t, 3> first = {{ 0, 0, 0 }};
    for (size_t i = 0, n = 1; i < size[0]; ++i) {
        for (size_t j = 0; j < size[1]; ++j) {
            for (size_t k = 0; k < size[2]; ++k, ++n) {
                array<size_t, 3> result = multi_range_find_if(
                    first
                  , size
                  , bind(&element_equal<double, 3>, cref(tensor), n, _1)
                );
                BOOST_CHECK_EQUAL(result[0], i);
                BOOST_CHECK_EQUAL(result[1], j);
                BOOST_CHECK_EQUAL(result[2], k);
            }
        }
    }
}

/**
 * find nonexistent value
 */
BOOST_FIXTURE_TEST_CASE( find_if_nonexistent, multi_array_fixture )
{
    array<size_t, 3> first = {{ 0, 0, 0 }};
    array<size_t, 3> result = multi_range_find_if(
        first
      , size
      , bind(&element_equal<double, 3>, cref(tensor), size[0] * size[1] * size[2] + 1, _1)
    );
    BOOST_CHECK_EQUAL(result[0], size[0]);
    BOOST_CHECK_EQUAL(result[1], size[1]);
    BOOST_CHECK_EQUAL(result[2], size[2]);
}


/**
 * find in empty range
 */
BOOST_FIXTURE_TEST_CASE( find_if_empty_range, multi_array_fixture )
{
    array<size_t, 3> first = {{ 0, 0, 0 }};
    array<size_t, 3> result = multi_range_find_if(
        first
      , first
      , bind(&element_equal<double, 3>, cref(tensor), size[0] * size[1] * size[2], _1)
    );
    BOOST_CHECK_EQUAL(result[0], first[0]);
    BOOST_CHECK_EQUAL(result[1], first[1]);
    BOOST_CHECK_EQUAL(result[2], first[2]);
}

template <typename T, size_t N>
void element_accumulate(multi_array<T, N> const& tensor, T& sum, array<size_t, N> const& index)
{
    sum += tensor(index);
}

/**
 * sum over all elements
 */
BOOST_FIXTURE_TEST_CASE( for_each_sum_elements, multi_array_fixture )
{
    double const N = size[0] * size[1] * size[2];
    double sum = 0;
    array<size_t, 3> first = {{ 0, 0, 0 }};
    array<size_t, 3> result = multi_range_for_each(
        first
      , size
      , bind(&element_accumulate<double, 3>, cref(tensor), ref(sum), _1)
    );
    BOOST_CHECK_EQUAL(result[0], size[0]);
    BOOST_CHECK_EQUAL(result[1], size[1]);
    BOOST_CHECK_EQUAL(result[2], size[2]);
    BOOST_CHECK_EQUAL(sum, N * (N + 1) / 2);
}

/**
 * sum over empty range
 */
BOOST_FIXTURE_TEST_CASE( for_each_sum_empty_range, multi_array_fixture )
{
    double sum = 0;
    array<size_t, 3> first = {{ 0, 0, 0 }};
    array<size_t, 3> result = multi_range_for_each(
        first
      , first
      , bind(&element_accumulate<double, 3>, cref(tensor), ref(sum), _1)
    );
    BOOST_CHECK_EQUAL(result[0], first[0]);
    BOOST_CHECK_EQUAL(result[1], first[1]);
    BOOST_CHECK_EQUAL(result[2], first[2]);
    BOOST_CHECK_EQUAL(sum, 0);
}
