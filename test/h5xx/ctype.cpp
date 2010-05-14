/*
 * Copyright © 2010  Peter Colberg
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ctype
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <vector>

#include <h5xx/ctype.hpp>

using namespace std;

typedef vector<pair<hid_t, hid_t> > type_pairs_vector;

static void add_native_types(type_pairs_vector& types)
{
    types.push_back(make_pair(h5xx::ctype<float>(), H5T_NATIVE_FLOAT));
    types.push_back(make_pair(h5xx::ctype<double>(), H5T_NATIVE_DOUBLE));
    types.push_back(make_pair(h5xx::ctype<signed char>(), H5T_NATIVE_CHAR));
    types.push_back(make_pair(h5xx::ctype<unsigned char>(), H5T_NATIVE_UCHAR));
    types.push_back(make_pair(h5xx::ctype<short>(), H5T_NATIVE_SHORT));
    types.push_back(make_pair(h5xx::ctype<unsigned short>(), H5T_NATIVE_USHORT));
    types.push_back(make_pair(h5xx::ctype<int>(), H5T_NATIVE_INT));
    types.push_back(make_pair(h5xx::ctype<unsigned int>(), H5T_NATIVE_UINT));
    types.push_back(make_pair(h5xx::ctype<long>(), H5T_NATIVE_LONG));
    types.push_back(make_pair(h5xx::ctype<unsigned long>(), H5T_NATIVE_ULONG));
}

static void add_integer_types(type_pairs_vector& types)
{
    types.push_back(make_pair(h5xx::ctype<int8_t>(), H5T_NATIVE_INT8));
    types.push_back(make_pair(h5xx::ctype<uint8_t>(), H5T_NATIVE_UINT8));
    types.push_back(make_pair(h5xx::ctype<int16_t>(), H5T_NATIVE_INT16));
    types.push_back(make_pair(h5xx::ctype<uint16_t>(), H5T_NATIVE_UINT16));
    types.push_back(make_pair(h5xx::ctype<int32_t>(), H5T_NATIVE_INT32));
    types.push_back(make_pair(h5xx::ctype<uint32_t>(), H5T_NATIVE_UINT32));
    types.push_back(make_pair(h5xx::ctype<int64_t>(), H5T_NATIVE_INT64));
    types.push_back(make_pair(h5xx::ctype<uint64_t>(), H5T_NATIVE_UINT64));
}

static void test_equality(type_pairs_vector const& types)
{
    for (size_t i = 0; i < types.size(); ++i) {
        BOOST_CHECK(0 < H5Tequal(types[i].first, types[i].second));
    }
}

static void test_inequality(type_pairs_vector const& types)
{
    for (size_t i = 0; i < types.size(); ++i) {
        for (size_t j = 0; j < types.size(); ++j) {
            if (i != j) {
                BOOST_CHECK(0 == H5Tequal(types[i].first, types[j].second));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_native_types_equality)
{
    type_pairs_vector types;
    add_native_types(types);
    test_equality(types);
}

BOOST_AUTO_TEST_CASE(test_native_types_inequality)
{
    type_pairs_vector types;
    add_native_types(types);
    test_inequality(types);
}

BOOST_AUTO_TEST_CASE(test_integer_types_equality)
{
    type_pairs_vector types;
    add_integer_types(types);
    test_equality(types);
}

BOOST_AUTO_TEST_CASE(test_integer_types_inequality)
{
    type_pairs_vector types;
    add_integer_types(types);
    test_inequality(types);
}
