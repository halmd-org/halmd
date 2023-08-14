/*
 * Copyright © 2011-2012  Peter Colberg
 * Copyright © 2020 Jaslo Ziska
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

#define BOOST_TEST_MODULE multi_range
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <boost/assign.hpp>
#include <boost/bind/bind.hpp>
#include <boost/multi_array.hpp>
#include <iterator> // std::back_inserter

#include <halmd/algorithm/multi_range.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace std;

template <typename T, size_t N>
bool element_equal(multi_array<T, N> const& tensor, T const& value, boost::array<size_t, N> const& index)
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

    boost::array<size_t, 3> const size;
    multi_array<double, 3> tensor;
};

/**
 * find each element by value
 */
BOOST_FIXTURE_TEST_CASE( find_if_every_element, multi_array_fixture )
{
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    for (size_t i = 0, n = 1; i < size[0]; ++i) {
        for (size_t j = 0; j < size[1]; ++j) {
            for (size_t k = 0; k < size[2]; ++k, ++n) {
                boost::array<size_t, 3> result = multi_range_find_if(
                    first
                  , size
                  , bind(&element_equal<double, 3>, boost::cref(tensor), n, boost::placeholders::_1)
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
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    boost::array<size_t, 3> result = multi_range_find_if(
        first
      , size
      , bind(&element_equal<double, 3>, boost::cref(tensor), size[0] * size[1] * size[2] + 1, boost::placeholders::_1)
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
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    boost::array<size_t, 3> result = multi_range_find_if(
        first
      , first
      , bind(&element_equal<double, 3>, boost::cref(tensor), size[0] * size[1] * size[2], boost::placeholders::_1)
    );
    BOOST_CHECK_EQUAL(result[0], first[0]);
    BOOST_CHECK_EQUAL(result[1], first[1]);
    BOOST_CHECK_EQUAL(result[2], first[2]);
}

template <typename T, size_t N>
class multi_array_accumulate
{
public:
    explicit multi_array_accumulate(multi_array<T, N> const& tensor) : tensor_(tensor), sum_(0) {}

    void operator()(boost::array<size_t, N> const& index)
    {
        sum_ += tensor_(index);
    }

    T operator()() const
    {
        return sum_;
    }

private:
    multi_array<T, N> const& tensor_;
    T sum_;
};

/**
 * sum over all elements
 */
BOOST_FIXTURE_TEST_CASE( for_each_sum_elements, multi_array_fixture )
{
    double const N = size[0] * size[1] * size[2];
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    multi_array_accumulate<double, 3> result = multi_range_for_each(
        first
      , size
      , multi_array_accumulate<double, 3>(tensor)
    );
    BOOST_CHECK_EQUAL(result(), N * (N + 1) / 2);
}

/**
 * sum over empty range
 */
BOOST_FIXTURE_TEST_CASE( for_each_sum_empty_range, multi_array_fixture )
{
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    multi_array_accumulate<double, 3> result = multi_range_for_each(
        first
      , first
      , multi_array_accumulate<double, 3>(tensor)
    );
    BOOST_CHECK_EQUAL(result(), 0);
}

template <typename multi_array_type, typename iterator_type>
class multi_range_copy
{
public:
    multi_range_copy(multi_array_type const& input, iterator_type output) : input_(input), output_(output) {}

    /**
     * Copy element from multi_array to output iterator.
     */
    template <typename index_type>
    bool operator()(index_type const& index)
    {
        *output_++ = input_(index);
        return false;
    }

private:
    multi_array_type const& input_;
    iterator_type output_;
};

template <typename multi_array_type, typename iterator_type>
static multi_range_copy<multi_array_type, iterator_type>
make_multi_range_copy(multi_array_type const& input, iterator_type output)
{
    return multi_range_copy<multi_array_type, iterator_type>(input, output);
}

/**
 * Check multi_range_for_each iterates in C storage order,
 * which is the default storage order of multi_array.
 */
BOOST_FIXTURE_TEST_CASE( multi_range_for_each_iteration_order , multi_array_fixture )
{
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    vector<double> result;
    result.reserve(tensor.num_elements());
    multi_range_for_each(
        first
      , size
      , make_multi_range_copy(tensor, back_inserter(result))
    );
    // check result vector has same order as multi_array storage
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result.begin()
      , result.end()
      , tensor.data()
      , tensor.data() + tensor.num_elements()
    );
}

/**
 * Check multi_range_find_if iterates in C storage order,
 * which is the default storage order of multi_array.
 */
BOOST_FIXTURE_TEST_CASE( multi_range_find_if_iteration_order , multi_array_fixture )
{
    boost::array<size_t, 3> first = {{ 0, 0, 0 }};
    vector<double> result;
    result.reserve(tensor.num_elements());
    multi_range_find_if(
        first
      , size
      , make_multi_range_copy(tensor, back_inserter(result))
    );
    // check result vector has same order as multi_array storage
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result.begin()
      , result.end()
      , tensor.data()
      , tensor.data() + tensor.num_elements()
    );
}

struct multi_array_4_fixture
{
    typedef multi_array<double, 4> array_type;
    typedef boost::array<size_t, 4> index_type;

    index_type const size;
    array_type tensor;

    /**
     * Initialise elements of multi_array with ascending sequence of integers.
     * The sequence is continuous in memory, following C storage order.
     */
    explicit multi_array_4_fixture(index_type const& size) : size(size), tensor(size)
    {
        size_t n = 0;
        for (size_t i = 0; i < size[0]; ++i) {
            for (size_t j = 0; j < size[1]; ++j) {
                for (size_t k = 0; k < size[2]; ++k) {
                    for (size_t l = 0; l < size[3]; ++l) {
                        tensor[i][j][k][l] = ++n;
                    }
                }
            }
        }
    }
};

/**
 * Compute sum of elements of multi_array.
 */
template <typename input_type>
struct multi_array_sum
{
    explicit multi_array_sum(input_type const& input) : input_(input), sum_(0) {}

    /**
     * Add multi_array element to sum.
     */
    template <typename index_type>
    void operator()(index_type const& index)
    {
        sum_ += input_(index);
    }

    /**
     * Returns sum.
     */
    double operator()() const
    {
        return sum_;
    }

private:
    input_type const& input_;
    double sum_;
};

/**
 * Measure average time of summing multi_array using multi_range_for_each.
 */
struct time_multi_range_for_each : multi_array_4_fixture
{
    explicit time_multi_range_for_each(index_type const& size) : multi_array_4_fixture(size) {}

    void operator()(size_t iterations) const
    {
        BOOST_TEST_MESSAGE( "  sum " << size[0] << "×" << size[1] << "×" << size[2] << "×" << size[3] << " elements using multi_range_for_each");
        BOOST_TEST_MESSAGE( "  warm up " << iterations << " iterations");
        accumulator<double> acc;
        for (size_t i = 0; i < iterations; ++i) {
            index_type first = {{ 0, 0, 0, 0 }};
            multi_array_sum<array_type> result = multi_range_for_each(
                first
              , size
              , multi_array_sum<array_type>(tensor)
            );
            acc(result());
        }

        BOOST_TEST_MESSAGE( "  measure " << iterations << " iterations" );
        accumulator<double> elapsed;
        for (size_t i = 0; i < iterations; ++i) {
            scoped_timer<timer> t(elapsed);
            index_type first = {{ 0, 0, 0, 0 }};
            multi_array_sum<array_type> result = multi_range_for_each(
                first
              , size
              , multi_array_sum<array_type>(tensor)
            );
            acc(result());
        }

        uint64_t const n = tensor.num_elements();
        BOOST_TEST_MESSAGE( "  elapsed time " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
        BOOST_CHECK_EQUAL( count(elapsed), iterations );
        BOOST_CHECK_EQUAL( mean(acc), n * (n + 1) / 2 );
        BOOST_CHECK_EQUAL( error_of_mean(acc), 0 );
        BOOST_CHECK_EQUAL( count(acc), 2 * iterations );
    }
};

/**
 * Measure average time of summing multi_array using nested for loops.
 */
struct time_nested_for_loops : multi_array_4_fixture
{
    explicit time_nested_for_loops(index_type const& size) : multi_array_4_fixture(size) {}

    void operator()(size_t iterations) const
    {
        BOOST_TEST_MESSAGE( "  sum " << size[0] << "×" << size[1] << "×" << size[2] << "×" << size[3] << " elements using nested for loops");
        BOOST_TEST_MESSAGE( "  warm up " << iterations << " iterations");
        accumulator<double> acc;
        for (size_t i = 0; i < iterations; ++i) {
            double sum = 0;
            index_type k;
            for (k[0] = 0; k[0] != size[0]; ++k[0]) {
                for (k[1] = 0; k[1] != size[1]; ++k[1]) {
                    for (k[2] = 0; k[2] != size[2]; ++k[2]) {
                        for (k[3] = 0; k[3] != size[3]; ++k[3]) {
                            sum += tensor(k);
                        }
                    }
                }
            }
            acc(sum);
        }

        BOOST_TEST_MESSAGE( "  measure " << iterations << " iterations" );
        accumulator<double> elapsed;
        for (size_t i = 0; i < iterations; ++i) {
            scoped_timer<timer> t(elapsed);
            double sum = 0;
            index_type k;
            for (k[0] = 0; k[0] != size[0]; ++k[0]) {
                for (k[1] = 0; k[1] != size[1]; ++k[1]) {
                    for (k[2] = 0; k[2] != size[2]; ++k[2]) {
                        for (k[3] = 0; k[3] != size[3]; ++k[3]) {
                            sum += tensor(k);
                        }
                    }
                }
            }
            acc(sum);
        }

        uint64_t const n = tensor.num_elements();
        BOOST_TEST_MESSAGE( "  elapsed time " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
        BOOST_CHECK_EQUAL( count(elapsed), iterations );
        BOOST_CHECK_EQUAL( mean(acc), n * (n + 1) / 2 );
        BOOST_CHECK_EQUAL( error_of_mean(acc), 0 );
        BOOST_CHECK_EQUAL( count(acc), 2 * iterations );
    }
};

/**
 * Data-driven test case regristration.
 */
using namespace boost::unit_test;

auto sizes_range = data::xrange(20, 101, 20);

BOOST_DATA_TEST_CASE( test_time_multi_range_for_each, sizes_range, size ) {
    time_multi_range_for_each::index_type sizes = list_of(size)(size)(size)(size);
    time_multi_range_for_each test(sizes);
    test(10);
}
BOOST_DATA_TEST_CASE( test_time_nested_for_loops, sizes_range, size ) {
    time_nested_for_loops::index_type sizes = list_of(size)(size)(size)(size);
    time_nested_for_loops test(sizes);
    test(10);
}
