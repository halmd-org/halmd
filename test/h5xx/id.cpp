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
#define BOOST_TEST_MODULE test_h5xx_id
#include <boost/test/unit_test.hpp>

// #define H5XX_DEBUG
#include <h5xx/id.hpp>

#define H5XX_TEST_FILENAME "test_h5xx_id.h5"

using namespace std;

/**
 * expose h5xx::id constructor
 *
 * Call this only *once* per identifier created with H5 C function.
 */
struct _expose_id : public h5xx::id
{
    _expose_id(hid_t id) : h5xx::id(id) {}
};

/**
 * returns identifier reference count
 */
inline int count(hid_t id)
{
    int count;
    H5XX_CHECK(count = H5Iget_ref(id));
    return count;
}

BOOST_AUTO_TEST_CASE(test_id)
{
    // construct two h5xx::id from file identifier
    hid_t file;
    H5XX_CHECK(file = H5Fcreate(H5XX_TEST_FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    BOOST_CHECK(count(file) == 1);
    h5xx::id id1 = _expose_id(file);
    BOOST_CHECK(count(file) == 1);
    h5xx::id id2 = id1;
    BOOST_CHECK(count(file) == 2);
    BOOST_CHECK(id1 == id2);
    BOOST_CHECK(!(id1 != id2));

    // construct two h5xx::id from dataspace identifier
    hid_t dataspace;
    H5XX_CHECK(dataspace = H5Screate(H5S_SCALAR));
    h5xx::id id3 = _expose_id(dataspace);
    BOOST_CHECK(count(dataspace) == 1);
    BOOST_CHECK(id1 != id3);
    BOOST_CHECK(!(id1 == id3));
    h5xx::id id4(id3);
    BOOST_CHECK(count(dataspace) == 2);
    BOOST_CHECK(id3 == id4);

    // mingle file and dataspace identifiers
    id1 = id3;
    BOOST_CHECK(count(file) == 1);
    BOOST_CHECK(count(dataspace) == 3);
    BOOST_CHECK(id1 == id3);
    BOOST_CHECK(id1 != id2);
    id4 = id2;
    BOOST_CHECK(count(file) == 2);
    BOOST_CHECK(count(dataspace) == 2);
    BOOST_CHECK(id4 == id2);
    BOOST_CHECK(id4 != id3);
    id2 = id4 = id1;
    BOOST_CHECK(count(dataspace) == 4);
    BOOST_CHECK(id2 == id4);
    BOOST_CHECK(id2 == id1);
    BOOST_CHECK(id2 == id3);
    try {
        count(file);
        BOOST_FAIL("failed to not get file identifier reference count");
    }
    catch (h5xx::error const& e) {
        BOOST_CHECK(e.count(make_pair(H5E_ATOM, H5E_BADATOM)));
    }
}
