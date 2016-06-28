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

#define BOOST_TEST_MODULE index
#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>

#include <halmd/utility/multi_index.hpp>
#include <test/tools/ctest.hpp>

using namespace halmd;

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

BOOST_AUTO_TEST_CASE( multi_index_1_dim )
{
    boost::array<size_t, 1> index = {{ 17 }};
    boost::array<size_t, 1> dims = {{ 19 }};
    size_t offset = 17;
    BOOST_CHECK_EQUAL( multi_index_to_offset(index, dims), offset );
    BOOST_CHECK_EQUAL( offset_to_multi_index(offset, dims), index );
}

BOOST_AUTO_TEST_CASE( multi_index_2_dim )
{
    boost::array<size_t, 2> index = {{ 17, 11 }};
    boost::array<size_t, 2> dims = {{ 19, 13 }};
    size_t offset = 17 + 11 * 19;
    BOOST_CHECK_EQUAL( multi_index_to_offset(index, dims), offset );
    BOOST_CHECK_EQUAL( offset_to_multi_index(offset, dims), index );
}

BOOST_AUTO_TEST_CASE( multi_index_3_dim )
{
    boost::array<size_t, 3> index = {{ 11, 17, 13 }};
    boost::array<size_t, 3> dims = {{ 19, 29, 23 }};
    size_t offset = 11 + 17 * 19 + 13 * (19 * 29);
    BOOST_CHECK_EQUAL( multi_index_to_offset(index, dims), offset );
    BOOST_CHECK_EQUAL( offset_to_multi_index(offset, dims), index );
}

BOOST_AUTO_TEST_CASE( multi_index_4_dim )
{
    boost::array<size_t, 4> index = {{ 11, 79, 17, 13 }};
    boost::array<size_t, 4> dims = {{ 19, 131, 29, 23 }};
    size_t offset = 11 + 79 * 19 + 17 * (19 * 131) + 13 * (19 * 131 * 29);
    BOOST_CHECK_EQUAL( multi_index_to_offset(index, dims), offset );
    BOOST_CHECK_EQUAL( offset_to_multi_index(offset, dims), index );
}

BOOST_AUTO_TEST_CASE( multi_index_5_dim )
{
    boost::array<size_t, 5> index = {{ 11, 79, 127, 17, 13 }};
    boost::array<size_t, 5> dims = {{ 19, 131, 137, 29, 23 }};
    size_t offset = 11 + 79 * 19 + 127 * (19 * 131) + 17 * (19 * 131 * 137) + 13 * (19 * 131 * 137 * 29);
    BOOST_CHECK_EQUAL( multi_index_to_offset(index, dims), offset );
    BOOST_CHECK_EQUAL( offset_to_multi_index(offset, dims), index );
}
