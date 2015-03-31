/*
 * Copyright © 2013,2015 Felix Höfling
 * Copyright © 2015      Nicolas Höft
 * Copyright © 2012      Peter Colberg
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
#include <test/tools/ctest.hpp>
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
        decltype(begin(*static_cast<array_type*>(0)))
      , typename array_type::iterator
    >::value, "begin(…) returns iterator");

    static_assert(std::is_same<
        decltype(begin(*static_cast<array_type const*>(0)))
      , typename array_type::const_iterator
    >::value, "begin(…) returns const iterator");

    static_assert(std::is_same<
        decltype(cbegin(*static_cast<array_type const*>(0)))
      , typename array_type::const_iterator
    >::value, "cbegin(…) returns const_iterator");

    static_assert(std::is_same<
        decltype((*static_cast<array_type*>(0)).end())
      , typename array_type::iterator
    >::value, "end() returns iterator");

    static_assert(std::is_same<
        decltype((*static_cast<array_type const*>(0)).end())
      , typename array_type::const_iterator
    >::value, "end() const returns const_iterator");

    static_assert(std::is_same<
        decltype(end(*static_cast<array_type*>(0)))
      , typename array_type::iterator
    >::value, "end(…) returns iterator");

    static_assert(std::is_same<
        decltype(end(*static_cast<array_type const*>(0)))
      , typename array_type::const_iterator
    >::value, "end(…) returns const iterator");

    static_assert(std::is_same<
        decltype(cend(*static_cast<array_type const*>(0)))
      , typename array_type::const_iterator
    >::value, "cend(…) returns const_iterator");

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
static void test_move_constructor(halmd::raw_array<T>& array)
{
    typedef halmd::raw_array<T> array_type;
    typedef typename array_type::size_type size_type;
    size_type const size = array.size();
    size_type const capacity = array.capacity();
    typename array_type::pointer const begin = &*array.begin();

    BOOST_CHECK_EQUAL( array.size(), size );
    BOOST_CHECK_EQUAL( array.capacity(), capacity );
    BOOST_CHECK_EQUAL( &*array.begin(), begin );
    BOOST_CHECK_EQUAL( &*array.end(), begin + size );

    // construct array from xvalue reference of non-empty array
    array_type array2(std::move(array));

    BOOST_CHECK_EQUAL( array.size(), 0lu );
    BOOST_CHECK_EQUAL( array.capacity(), 0lu );
    BOOST_CHECK_EQUAL( array.begin(), array.end() );
    BOOST_CHECK_EQUAL( array2.size(), size );
    BOOST_CHECK_EQUAL( array2.capacity(), capacity );
    BOOST_CHECK_EQUAL( &*array2.begin(), begin );
    BOOST_CHECK_EQUAL( &*array2.end(), begin + size );

    // construct array from xvalue reference of empty array
    array_type array3(std::move(array));

    BOOST_CHECK_EQUAL( array.size(), 0lu );
    BOOST_CHECK_EQUAL( array.capacity(), 0lu );
    BOOST_CHECK_EQUAL( array.begin(), array.end() );
    BOOST_CHECK_EQUAL( array3.size(), 0lu );
    BOOST_CHECK_EQUAL( array3.capacity(), 0lu );
    BOOST_CHECK_EQUAL( array3.begin(), array3.end() );

    // array_type array4(array);   // must not compile due to deleted copy constructor
}

template <typename T>
static void test_move_assignment(halmd::raw_array<T>& array)
{
    typedef halmd::raw_array<T> array_type;
    typedef typename array_type::size_type size_type;
    size_type const size = array.size();
    size_type const capacity = array.capacity();
    typename array_type::pointer const begin = &*array.begin();

    BOOST_CHECK_EQUAL( array.size(), size );
    BOOST_CHECK_EQUAL( array.capacity(), capacity );
    BOOST_CHECK_EQUAL( &*array.begin(), begin );
    BOOST_CHECK_EQUAL( &*array.end(), begin + size );

    // construct non-empty array
    size_type size2 = size / 2;
    size_type capacity2 = size;
    array_type array2(size2);
    array2.reserve(capacity2);
    typename array_type::pointer const begin2 = &*array2.begin();
    BOOST_CHECK_EQUAL( array2.size(), size2 );
    BOOST_CHECK_EQUAL( array2.capacity(), capacity2 );
    BOOST_CHECK_EQUAL( array2.end(), begin2 + size2);

    // swap non-empty arrays
    swap(array, array2);
    BOOST_CHECK_EQUAL( array.size(), size2 );
    BOOST_CHECK_EQUAL( array.capacity(), capacity2 );
    BOOST_CHECK_EQUAL( array.begin(), begin2);
    BOOST_CHECK_EQUAL( array.end(), begin2 + size2 );
    BOOST_CHECK_EQUAL( array2.size(), size );
    BOOST_CHECK_EQUAL( array2.capacity(), size );
    BOOST_CHECK_EQUAL( array2.begin(), begin);
    BOOST_CHECK_EQUAL( array2.end(), begin + size );

    // assign xvalue reference of non-empty array
    array2 = std::move(array);
    // array2 = array;   // must not compile due to deleted copy constructor

    BOOST_CHECK_EQUAL( array.size(), 0lu );
    BOOST_CHECK_EQUAL( array.capacity(), 0lu );
    BOOST_CHECK_EQUAL( array.begin(), array.end() );
    BOOST_CHECK_EQUAL( array2.size(), size2 );
    BOOST_CHECK_EQUAL( array2.capacity(), capacity2 );
    BOOST_CHECK_EQUAL( &*array2.begin(), begin2 );
    BOOST_CHECK_EQUAL( &*array2.end(), begin2 + size2 );

    // construct empty array and assign xvalue reference of empty array
    array_type array3;
    array = std::move(array3);

    BOOST_CHECK_EQUAL( array.size(), 0lu );
    BOOST_CHECK_EQUAL( array.capacity(), 0lu );
    BOOST_CHECK_EQUAL( array.begin(), array.end() );
    BOOST_CHECK_EQUAL( array3.size(), 0lu );
    BOOST_CHECK_EQUAL( array3.capacity(), 0lu );
    BOOST_CHECK_EQUAL( array3.begin(), array3.end() );

    // assign xvalue reference of self
    array2 = std::move(array2);

    BOOST_CHECK_EQUAL( array2.size(), size2 );
    BOOST_CHECK_EQUAL( array2.capacity(), capacity2 );
    BOOST_CHECK_EQUAL( &*array2.begin(), begin2 );
    BOOST_CHECK_EQUAL( &*array2.end(), begin2 + size2 );
}

template <typename T>
static void test_resize(halmd::raw_array<T>& array)
{
    typedef halmd::raw_array<T> array_type;
    typedef typename array_type::size_type size_type;

    size_type const size = array.size();
    size_type const capacity = array.capacity();

    std::copy(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
      , array.begin()
    );
    array.resize(size / 2);
    BOOST_CHECK_EQUAL( array.capacity(), capacity );
    BOOST_CHECK_EQUAL( array.size(), size / 2 );

    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(size + 1)
      , array.begin()
      , array.begin() + size
    );

    array.resize(size); // restore old size

    BOOST_CHECK_EQUAL( array.capacity(), capacity );
    BOOST_CHECK_EQUAL( array.size(), size );

    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(size + 1)
      , array.begin()
      , array.end()
    );

    // double the array size, old elements must be copied
    array.resize(size * 2);

    BOOST_CHECK_EQUAL( array.size(), size * 2 );
    BOOST_CHECK_GE( array.capacity(), size * 2 );

    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(size + 1)
      , array.begin()
      , array.begin() + size
    );
}

template <typename T>
static void test_reserve(halmd::raw_array<T>& array)
{
    typedef halmd::raw_array<T> array_type;
    typedef typename array_type::size_type size_type;

    size_type const size = array.size();
    size_type const capacity = array.capacity();

    std::copy(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(array.size() + 1)
      , array.begin()
    );

    // increase capacity
    size_type new_capacity = capacity * 2;
    array.reserve(new_capacity);
    BOOST_CHECK_EQUAL( array.size(), size );
    BOOST_CHECK_EQUAL( array.capacity(), new_capacity );

    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(size + 1)
      , array.begin()
      , array.end()
    );

    // reserve less memory than size(), should do nothing
    array.reserve(size / 2);
    BOOST_CHECK_EQUAL( array.size(), size );
    BOOST_CHECK_EQUAL( array.capacity(), new_capacity );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        boost::counting_iterator<T>(1)
      , boost::counting_iterator<T>(size + 1)
      , array.begin()
      , array.end()
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

    auto move_constructor = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_move_constructor(array);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( move_constructor ));

    auto move_assignment = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_move_assignment(array);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( move_assignment ));

    auto resize = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_resize(array);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( resize ));

    auto reserve = [=]() {
        BOOST_TEST_MESSAGE( " " << halmd::demangled_name<array_type>() << " of size " << size );
        array_type array(size);
        test_reserve(array);
    };
    framework::master_test_suite().add(BOOST_TEST_CASE( reserve ));
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
