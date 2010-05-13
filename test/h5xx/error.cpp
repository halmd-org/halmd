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
#define BOOST_TEST_MODULE test_h5xx_error
#include <boost/test/unit_test.hpp>
#include <cstring>

#include <h5xx/error.hpp>

BOOST_AUTO_TEST_CASE(test_success)
{
    hid_t ret, id;
    H5XX_CALL(ret = id = H5Screate(H5S_SCALAR));
    BOOST_REQUIRE(ret >= 0);
    H5XX_CALL(ret = H5Sclose(id));
    BOOST_REQUIRE(ret >= 0);
}

BOOST_AUTO_TEST_CASE(test_throw_exception)
{
    hid_t ret, id;
    H5XX_CALL(ret = id = H5Screate(H5S_SCALAR));
    H5XX_CALL(ret = H5Sclose(id));
    try {
        H5XX_CALL(ret = H5Sclose(id));
    }
    catch (h5xx::error const&) {
        BOOST_REQUIRE(ret < 0);
    }
}

BOOST_AUTO_TEST_CASE(test_disable_handler)
{
    H5E_auto_t auto_ = reinterpret_cast<H5E_auto_t>(H5Eprint);
    void* data_ = stderr;
    H5XX_CALL(H5Eget_auto(H5E_DEFAULT, &auto_, &data_));
    BOOST_REQUIRE(NULL == auto_);
    BOOST_REQUIRE(NULL == data_);
}

BOOST_AUTO_TEST_CASE(test_reset_default_error_handler)
{
    H5E_auto_t auto_ = NULL;
    void* data_ = stderr;
    H5Eget_auto(H5E_DEFAULT, &auto_, &data_);
    BOOST_REQUIRE(reinterpret_cast<H5E_auto_t>(H5Eprint) == auto_);
    BOOST_REQUIRE(NULL == data_);
}

BOOST_AUTO_TEST_CASE(test_library_error_description)
{
    hid_t ret, id;
    H5XX_CALL(ret = id = H5Screate(H5S_SCALAR));
    H5XX_CALL(ret = H5Sclose(id));
    try {
        H5XX_CALL(ret = H5Sclose(id));
    }
    catch (h5xx::error const& e) {
        BOOST_REQUIRE(0 == strcmp(e.what(), "H5Sclose: not a dataspace"));
    }
}

BOOST_AUTO_TEST_CASE(test_custom_error_description)
{
    try {
        throw h5xx::error("test custom error description");
    }
    catch (h5xx::error const& e) {
        BOOST_REQUIRE(0 == strcmp(e.what(), "test custom error description"));
    }
}
