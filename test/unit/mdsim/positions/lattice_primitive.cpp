/*
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

#define BOOST_TEST_MODULE lattice_primitive
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>

struct hexagonal_lattice_fixture
{
    typedef halmd::close_packed_lattice<boost::array<float, 2>, boost::array<size_t, 2> > lattice_type;
    typedef lattice_type::result_type result_type;
    typedef lattice_type::shape_type shape_type;

    static shape_type const& shape()
    {
        static shape_type shape = {{4, 3}};
        return shape;
    }

    static boost::array<result_type, 24> const& result()
    {
        static boost::array<result_type, 24> result = {{
            {{0.25, 0.25}}, {{0.75, 0.75}}
          , {{1.25, 0.25}}, {{1.75, 0.75}}
          , {{2.25, 0.25}}, {{2.75, 0.75}}
          , {{3.25, 0.25}}, {{3.75, 0.75}}
          , {{0.25, 1.25}}, {{0.75, 1.75}}
          , {{1.25, 1.25}}, {{1.75, 1.75}}
          , {{2.25, 1.25}}, {{2.75, 1.75}}
          , {{3.25, 1.25}}, {{3.75, 1.75}}
          , {{0.25, 2.25}}, {{0.75, 2.75}}
          , {{1.25, 2.25}}, {{1.75, 2.75}}
          , {{2.25, 2.25}}, {{2.75, 2.75}}
          , {{3.25, 2.25}}, {{3.75, 2.75}}
        }};
        return result;
    }
};

struct face_centered_cubic_lattice_fixture
{
    typedef halmd::close_packed_lattice<boost::array<float, 3>, boost::array<size_t, 3> > lattice_type;
    typedef lattice_type::result_type result_type;
    typedef lattice_type::shape_type shape_type;

    static shape_type const& shape()
    {
        static shape_type shape = {{4, 2, 3}};
        return shape;
    }

    static boost::array<result_type, 96> const& result()
    {
        static boost::array<result_type, 96> result = {{
            {{0.25, 0.25, 0.25}}, {{0.75, 0.75, 0.25}}, {{0.75, 0.25, 0.75}}, {{0.25, 0.75, 0.75}}
          , {{1.25, 0.25, 0.25}}, {{1.75, 0.75, 0.25}}, {{1.75, 0.25, 0.75}}, {{1.25, 0.75, 0.75}}
          , {{2.25, 0.25, 0.25}}, {{2.75, 0.75, 0.25}}, {{2.75, 0.25, 0.75}}, {{2.25, 0.75, 0.75}}
          , {{3.25, 0.25, 0.25}}, {{3.75, 0.75, 0.25}}, {{3.75, 0.25, 0.75}}, {{3.25, 0.75, 0.75}}
          , {{0.25, 1.25, 0.25}}, {{0.75, 1.75, 0.25}}, {{0.75, 1.25, 0.75}}, {{0.25, 1.75, 0.75}}
          , {{1.25, 1.25, 0.25}}, {{1.75, 1.75, 0.25}}, {{1.75, 1.25, 0.75}}, {{1.25, 1.75, 0.75}}
          , {{2.25, 1.25, 0.25}}, {{2.75, 1.75, 0.25}}, {{2.75, 1.25, 0.75}}, {{2.25, 1.75, 0.75}}
          , {{3.25, 1.25, 0.25}}, {{3.75, 1.75, 0.25}}, {{3.75, 1.25, 0.75}}, {{3.25, 1.75, 0.75}}
          , {{0.25, 0.25, 1.25}}, {{0.75, 0.75, 1.25}}, {{0.75, 0.25, 1.75}}, {{0.25, 0.75, 1.75}}
          , {{1.25, 0.25, 1.25}}, {{1.75, 0.75, 1.25}}, {{1.75, 0.25, 1.75}}, {{1.25, 0.75, 1.75}}
          , {{2.25, 0.25, 1.25}}, {{2.75, 0.75, 1.25}}, {{2.75, 0.25, 1.75}}, {{2.25, 0.75, 1.75}}
          , {{3.25, 0.25, 1.25}}, {{3.75, 0.75, 1.25}}, {{3.75, 0.25, 1.75}}, {{3.25, 0.75, 1.75}}
          , {{0.25, 1.25, 1.25}}, {{0.75, 1.75, 1.25}}, {{0.75, 1.25, 1.75}}, {{0.25, 1.75, 1.75}}
          , {{1.25, 1.25, 1.25}}, {{1.75, 1.75, 1.25}}, {{1.75, 1.25, 1.75}}, {{1.25, 1.75, 1.75}}
          , {{2.25, 1.25, 1.25}}, {{2.75, 1.75, 1.25}}, {{2.75, 1.25, 1.75}}, {{2.25, 1.75, 1.75}}
          , {{3.25, 1.25, 1.25}}, {{3.75, 1.75, 1.25}}, {{3.75, 1.25, 1.75}}, {{3.25, 1.75, 1.75}}
          , {{0.25, 0.25, 2.25}}, {{0.75, 0.75, 2.25}}, {{0.75, 0.25, 2.75}}, {{0.25, 0.75, 2.75}}
          , {{1.25, 0.25, 2.25}}, {{1.75, 0.75, 2.25}}, {{1.75, 0.25, 2.75}}, {{1.25, 0.75, 2.75}}
          , {{2.25, 0.25, 2.25}}, {{2.75, 0.75, 2.25}}, {{2.75, 0.25, 2.75}}, {{2.25, 0.75, 2.75}}
          , {{3.25, 0.25, 2.25}}, {{3.75, 0.75, 2.25}}, {{3.75, 0.25, 2.75}}, {{3.25, 0.75, 2.75}}
          , {{0.25, 1.25, 2.25}}, {{0.75, 1.75, 2.25}}, {{0.75, 1.25, 2.75}}, {{0.25, 1.75, 2.75}}
          , {{1.25, 1.25, 2.25}}, {{1.75, 1.75, 2.25}}, {{1.75, 1.25, 2.75}}, {{1.25, 1.75, 2.75}}
          , {{2.25, 1.25, 2.25}}, {{2.75, 1.75, 2.25}}, {{2.75, 1.25, 2.75}}, {{2.25, 1.75, 2.75}}
          , {{3.25, 1.25, 2.25}}, {{3.75, 1.75, 2.25}}, {{3.75, 1.25, 2.75}}, {{3.25, 1.75, 2.75}}
        }};
        return result;
    }
};

struct square_lattice_fixture
{
    typedef halmd::primitive_lattice<boost::array<float, 2>, boost::array<size_t, 2> > lattice_type;
    typedef lattice_type::result_type result_type;
    typedef lattice_type::shape_type shape_type;

    static shape_type const& shape()
    {
        static shape_type shape = {{4, 3}};
        return shape;
    }

    static boost::array<result_type, 12> const& result()
    {
        static boost::array<result_type, 12> result = {{
            {{0.5, 0.5}}, {{1.5, 0.5}}, {{2.5, 0.5}}, {{3.5, 0.5}}
          , {{0.5, 1.5}}, {{1.5, 1.5}}, {{2.5, 1.5}}, {{3.5, 1.5}}
          , {{0.5, 2.5}}, {{1.5, 2.5}}, {{2.5, 2.5}}, {{3.5, 2.5}}
        }};
        return result;
    }
};

struct primitive_cubic_lattice_fixture
{
    typedef halmd::primitive_lattice<boost::array<float, 3>, boost::array<size_t, 3> > lattice_type;
    typedef lattice_type::result_type result_type;
    typedef lattice_type::shape_type shape_type;

    static shape_type const& shape()
    {
        static shape_type shape = {{4, 2, 3}};
        return shape;
    }

    static boost::array<result_type, 24> const& result()
    {
        static boost::array<result_type, 24> result = {{
            {{0.5, 0.5, 0.5}}, {{1.5, 0.5, 0.5}}, {{2.5, 0.5, 0.5}}, {{3.5, 0.5, 0.5}}
          , {{0.5, 1.5, 0.5}}, {{1.5, 1.5, 0.5}}, {{2.5, 1.5, 0.5}}, {{3.5, 1.5, 0.5}}
          , {{0.5, 0.5, 1.5}}, {{1.5, 0.5, 1.5}}, {{2.5, 0.5, 1.5}}, {{3.5, 0.5, 1.5}}
          , {{0.5, 1.5, 1.5}}, {{1.5, 1.5, 1.5}}, {{2.5, 1.5, 1.5}}, {{3.5, 1.5, 1.5}}
          , {{0.5, 0.5, 2.5}}, {{1.5, 0.5, 2.5}}, {{2.5, 0.5, 2.5}}, {{3.5, 0.5, 2.5}}
          , {{0.5, 1.5, 2.5}}, {{1.5, 1.5, 2.5}}, {{2.5, 1.5, 2.5}}, {{3.5, 1.5, 2.5}}
        }};
        return result;
    }
};

struct primitive_tesseractic_lattice_fixture
{
    typedef halmd::primitive_lattice<boost::array<float, 4>, boost::array<size_t, 4> > lattice_type;
    typedef lattice_type::result_type result_type;
    typedef lattice_type::shape_type shape_type;

    static shape_type const& shape()
    {
        static shape_type shape = {{4, 2, 3, 5}};
        return shape;
    }

    static boost::array<result_type, 120> const& result()
    {
        static boost::array<result_type, 120> result = {{
            {{0.5, 0.5, 0.5, 0.5}}, {{1.5, 0.5, 0.5, 0.5}}, {{2.5, 0.5, 0.5, 0.5}}, {{3.5, 0.5, 0.5, 0.5}}
          , {{0.5, 1.5, 0.5, 0.5}}, {{1.5, 1.5, 0.5, 0.5}}, {{2.5, 1.5, 0.5, 0.5}}, {{3.5, 1.5, 0.5, 0.5}}
          , {{0.5, 0.5, 1.5, 0.5}}, {{1.5, 0.5, 1.5, 0.5}}, {{2.5, 0.5, 1.5, 0.5}}, {{3.5, 0.5, 1.5, 0.5}}
          , {{0.5, 1.5, 1.5, 0.5}}, {{1.5, 1.5, 1.5, 0.5}}, {{2.5, 1.5, 1.5, 0.5}}, {{3.5, 1.5, 1.5, 0.5}}
          , {{0.5, 0.5, 2.5, 0.5}}, {{1.5, 0.5, 2.5, 0.5}}, {{2.5, 0.5, 2.5, 0.5}}, {{3.5, 0.5, 2.5, 0.5}}
          , {{0.5, 1.5, 2.5, 0.5}}, {{1.5, 1.5, 2.5, 0.5}}, {{2.5, 1.5, 2.5, 0.5}}, {{3.5, 1.5, 2.5, 0.5}}
          , {{0.5, 0.5, 0.5, 1.5}}, {{1.5, 0.5, 0.5, 1.5}}, {{2.5, 0.5, 0.5, 1.5}}, {{3.5, 0.5, 0.5, 1.5}}
          , {{0.5, 1.5, 0.5, 1.5}}, {{1.5, 1.5, 0.5, 1.5}}, {{2.5, 1.5, 0.5, 1.5}}, {{3.5, 1.5, 0.5, 1.5}}
          , {{0.5, 0.5, 1.5, 1.5}}, {{1.5, 0.5, 1.5, 1.5}}, {{2.5, 0.5, 1.5, 1.5}}, {{3.5, 0.5, 1.5, 1.5}}
          , {{0.5, 1.5, 1.5, 1.5}}, {{1.5, 1.5, 1.5, 1.5}}, {{2.5, 1.5, 1.5, 1.5}}, {{3.5, 1.5, 1.5, 1.5}}
          , {{0.5, 0.5, 2.5, 1.5}}, {{1.5, 0.5, 2.5, 1.5}}, {{2.5, 0.5, 2.5, 1.5}}, {{3.5, 0.5, 2.5, 1.5}}
          , {{0.5, 1.5, 2.5, 1.5}}, {{1.5, 1.5, 2.5, 1.5}}, {{2.5, 1.5, 2.5, 1.5}}, {{3.5, 1.5, 2.5, 1.5}}
          , {{0.5, 0.5, 0.5, 2.5}}, {{1.5, 0.5, 0.5, 2.5}}, {{2.5, 0.5, 0.5, 2.5}}, {{3.5, 0.5, 0.5, 2.5}}
          , {{0.5, 1.5, 0.5, 2.5}}, {{1.5, 1.5, 0.5, 2.5}}, {{2.5, 1.5, 0.5, 2.5}}, {{3.5, 1.5, 0.5, 2.5}}
          , {{0.5, 0.5, 1.5, 2.5}}, {{1.5, 0.5, 1.5, 2.5}}, {{2.5, 0.5, 1.5, 2.5}}, {{3.5, 0.5, 1.5, 2.5}}
          , {{0.5, 1.5, 1.5, 2.5}}, {{1.5, 1.5, 1.5, 2.5}}, {{2.5, 1.5, 1.5, 2.5}}, {{3.5, 1.5, 1.5, 2.5}}
          , {{0.5, 0.5, 2.5, 2.5}}, {{1.5, 0.5, 2.5, 2.5}}, {{2.5, 0.5, 2.5, 2.5}}, {{3.5, 0.5, 2.5, 2.5}}
          , {{0.5, 1.5, 2.5, 2.5}}, {{1.5, 1.5, 2.5, 2.5}}, {{2.5, 1.5, 2.5, 2.5}}, {{3.5, 1.5, 2.5, 2.5}}
          , {{0.5, 0.5, 0.5, 3.5}}, {{1.5, 0.5, 0.5, 3.5}}, {{2.5, 0.5, 0.5, 3.5}}, {{3.5, 0.5, 0.5, 3.5}}
          , {{0.5, 1.5, 0.5, 3.5}}, {{1.5, 1.5, 0.5, 3.5}}, {{2.5, 1.5, 0.5, 3.5}}, {{3.5, 1.5, 0.5, 3.5}}
          , {{0.5, 0.5, 1.5, 3.5}}, {{1.5, 0.5, 1.5, 3.5}}, {{2.5, 0.5, 1.5, 3.5}}, {{3.5, 0.5, 1.5, 3.5}}
          , {{0.5, 1.5, 1.5, 3.5}}, {{1.5, 1.5, 1.5, 3.5}}, {{2.5, 1.5, 1.5, 3.5}}, {{3.5, 1.5, 1.5, 3.5}}
          , {{0.5, 0.5, 2.5, 3.5}}, {{1.5, 0.5, 2.5, 3.5}}, {{2.5, 0.5, 2.5, 3.5}}, {{3.5, 0.5, 2.5, 3.5}}
          , {{0.5, 1.5, 2.5, 3.5}}, {{1.5, 1.5, 2.5, 3.5}}, {{2.5, 1.5, 2.5, 3.5}}, {{3.5, 1.5, 2.5, 3.5}}
          , {{0.5, 0.5, 0.5, 4.5}}, {{1.5, 0.5, 0.5, 4.5}}, {{2.5, 0.5, 0.5, 4.5}}, {{3.5, 0.5, 0.5, 4.5}}
          , {{0.5, 1.5, 0.5, 4.5}}, {{1.5, 1.5, 0.5, 4.5}}, {{2.5, 1.5, 0.5, 4.5}}, {{3.5, 1.5, 0.5, 4.5}}
          , {{0.5, 0.5, 1.5, 4.5}}, {{1.5, 0.5, 1.5, 4.5}}, {{2.5, 0.5, 1.5, 4.5}}, {{3.5, 0.5, 1.5, 4.5}}
          , {{0.5, 1.5, 1.5, 4.5}}, {{1.5, 1.5, 1.5, 4.5}}, {{2.5, 1.5, 1.5, 4.5}}, {{3.5, 1.5, 1.5, 4.5}}
          , {{0.5, 0.5, 2.5, 4.5}}, {{1.5, 0.5, 2.5, 4.5}}, {{2.5, 0.5, 2.5, 4.5}}, {{3.5, 0.5, 2.5, 4.5}}
          , {{0.5, 1.5, 2.5, 4.5}}, {{1.5, 1.5, 2.5, 4.5}}, {{2.5, 1.5, 2.5, 4.5}}, {{3.5, 1.5, 2.5, 4.5}}
        }};
        return result;
    }
};

namespace std {

/**
 * Output boost::array to stream for BOOST_CHECK_EQUAL.
 */
template <typename T, size_t N>
ostream& operator<<(ostream& os, boost::array<T, N> const& index)
{
    for (size_t i = 0; i + 1 < N; ++i) {
        os << index[i] << " ";
    }
    return os << index[N - 1];
}

} // namespace std

template <typename Fixture>
struct test_lattice
{
    typedef typename Fixture::lattice_type lattice_type;

    /**
     * Test lattice primitive typedefs.
     */
    typedef typename lattice_type::result_type result_type;
    typedef typename lattice_type::shape_type shape_type;
    typedef typename lattice_type::size_type size_type;

    /**
     * Compare output of lattice primitive to lattice fixture.
     *
     * The lattice functor is constructed twice, once for the
     * begin iterator, and once for the end iterator, to test
     * that a call of the functor does not have side effects.
     *
     * This method also tests the .shape() and .size() methods.
     */
    void operator()() const
    {
        lattice_type const lattice(Fixture::shape());
        BOOST_CHECK_EQUAL_COLLECTIONS(
            boost::make_transform_iterator(boost::make_counting_iterator(size_type(0)), lattice)
          , boost::make_transform_iterator(boost::make_counting_iterator(lattice.size()), lattice_type(lattice.shape()))
          , Fixture::result().begin()
          , Fixture::result().end()
        );
    }
};

HALMD_TEST_INIT( lattice_primitive )
{
    using namespace boost::unit_test;

    auto hexagonal_lattice = test_lattice<hexagonal_lattice_fixture>();
    framework::master_test_suite().add(BOOST_TEST_CASE( hexagonal_lattice ));

    auto face_centered_cubic_lattice = test_lattice<face_centered_cubic_lattice_fixture>();
    framework::master_test_suite().add(BOOST_TEST_CASE( face_centered_cubic_lattice ));

    auto square_lattice = test_lattice<square_lattice_fixture>();
    framework::master_test_suite().add(BOOST_TEST_CASE( square_lattice ));

    auto primitive_cubic_lattice = test_lattice<primitive_cubic_lattice_fixture>();
    framework::master_test_suite().add(BOOST_TEST_CASE( primitive_cubic_lattice ));

    auto primitive_tesseractic_lattice = test_lattice<primitive_tesseractic_lattice_fixture>();
    framework::master_test_suite().add(BOOST_TEST_CASE( primitive_tesseractic_lattice ));
}
