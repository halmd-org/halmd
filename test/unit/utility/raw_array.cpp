/*
 * Copyright © 2012 Peter Colberg
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE raw_array
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <type_traits>

#include <halmd/utility/demangle.hpp>
#include <halmd/utility/raw_array.hpp>
#include <test/tools/init.hpp>

template <typename T>
static void test_construct(halmd::raw_array<T> const& array, std::size_t size)
{
    typedef halmd::raw_array<T> array_type;

    BOOST_CHECK_EQUAL( array.size(), size );
    BOOST_CHECK_EQUAL( array.end() - array.begin(), std::ptrdiff_t(size) );

    static_assert(std::is_same<
        typename array_type::size_type
      , std::size_t
    >::value, "size_type is std::size_t");

    static_assert(std::is_same<
        typename array_type::difference_type
      , std::ptrdiff_t
    >::value, "difference_type is std::ptrdiff_t");

    static_assert(std::is_same<
        typename array_type::value_type
      , T
    >::value, "value_type is T");

    static_assert(std::is_same<
        typename array_type::pointer
      , T*
    >::value, "pointer is T*");

    static_assert(std::is_same<
        typename array_type::const_pointer
      , T const*
    >::value, "const_pointer is T const*");

    static_assert(std::is_same<
        typename array_type::iterator
      , T*
    >::value, "iterator is T*");

    static_assert(std::is_same<
        typename array_type::const_iterator
      , T const*
    >::value, "const_iterator is T const*");

    static_assert(std::is_same<
        typename array_type::reference
      , T&
    >::value, "reference is T&");

    static_assert(std::is_same<
        typename array_type::const_reference
      , T const&
    >::value, "const_reference is T const&");

    static_assert(std::is_same<
        decltype((static_cast<array_type*>(0))->begin())
      , typename array_type::iterator
    >::value, "begin() returns iterator");

    static_assert(std::is_same<
        decltype((*static_cast<array_type const*>(0)).begin())
      , typename array_type::const_iterator
    >::value, "begin() const returns const_iterator");

    static_assert(std::is_same<
        decltype((*static_cast<array_type*>(0)).end())
      , typename array_type::iterator
    >::value, "end() returns iterator");

    static_assert(std::is_same<
        decltype((*static_cast<array_type const*>(0)).end())
      , typename array_type::const_iterator
    >::value, "end() const returns const_iterator");

    static_assert(std::is_same<
        decltype((*static_cast<array_type*>(0))[0])
      , typename array_type::reference
    >::value, "operator[] returns reference");

    static_assert(std::is_same<
        decltype((*static_cast<array_type const*>(0))[0])
      , typename array_type::const_reference
    >::value, "operator[] const returns const_reference");
}

template <typename T>
static void test_iterator(halmd::raw_array<T>& array)
{
    typedef halmd::raw_array<T> array_type;
    array_type const& const_array = array;

    std::copy(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
      , array.begin()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
      , const_array.begin()
      , const_array.end()
    );

    std::vector<T> const sequence(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
    );
    std::reverse(
        array.begin()
      , array.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        sequence.rbegin()
      , sequence.rend()
      , const_array.begin()
      , const_array.end()
    );
}

template <typename T>
static void test_array_subscript(halmd::raw_array<T>& array)
{
    typedef halmd::raw_array<T> array_type;
    array_type const& const_array = array;

    for (std::size_t i = 0; i < array.size(); ++i) {
        array[i] = i + 1;
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
      , const_array.begin()
      , const_array.end()
    );

    std::vector<T> const sequence(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
    );
    std::reverse(
        array.begin()
      , array.end()
    );
    auto const_array_at_index = [&](std::size_t i) {
        return const_array[i];
    };
    BOOST_CHECK_EQUAL_COLLECTIONS(
        sequence.rbegin()
      , sequence.rend()
      , boost::make_transform_iterator(boost::counting_iterator<std::size_t>(0), const_array_at_index)
      , boost::make_transform_iterator(boost::counting_iterator<std::size_t>(array.size()), const_array_at_index)
    );
}

template <typename T>
static void test_suite(std::size_t size)
{
    typedef halmd::raw_array<T> array_type;

    using namespace boost::unit_test;

    auto construct = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_construct(array, size);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( construct ));

    auto iterator = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_iterator(array);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( iterator ));

    auto array_subscript = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_array_subscript(array);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( array_subscript ));
}

/**
 * Manual test case registration.
 */
HALMD_TEST_INIT( raw_array )
{
    for (std::size_t size : {0, 1, 2, 5, 11, 23, 47, 191, 383, 6143, 786431}) {
        test_suite<unsigned int>(size);
        test_suite<float>(size);
        test_suite<double>(size);
    }
}
