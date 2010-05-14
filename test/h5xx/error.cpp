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
        BOOST_FAIL("no h5xx::error exception thrown");
    }
    catch (h5xx::error const&) {
        BOOST_REQUIRE(ret < 0);
    }
}

BOOST_AUTO_TEST_CASE(test_disable_handler)
{
    H5E_auto_t efunc = reinterpret_cast<H5E_auto_t>(H5Eprint);
    void* edata = stderr;
    H5XX_CALL(H5Eget_auto(&efunc, &edata));
    BOOST_REQUIRE(NULL == efunc);
    BOOST_REQUIRE(NULL == edata);
}

BOOST_AUTO_TEST_CASE(test_reset_default_error_handler)
{
    H5E_auto_t efunc = NULL;
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
    void* edata = NULL;
#else
    void* edata = stderr;
#endif
    H5Eget_auto(&efunc, &edata);
#ifdef H5_USE_16_API_DEFAULT
    BOOST_REQUIRE(reinterpret_cast<H5E_auto_t>(H5Eprint1) == efunc);
#else
    BOOST_REQUIRE(reinterpret_cast<H5E_auto_t>(H5Eprint) == efunc);
#endif
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
    BOOST_REQUIRE(stderr == edata);
#else
    BOOST_REQUIRE(NULL == edata);
#endif
}

BOOST_AUTO_TEST_CASE(test_library_error_description)
{
    hid_t ret, id;
    H5XX_CALL(ret = id = H5Screate(H5S_SCALAR));
    H5XX_CALL(ret = H5Sclose(id));
    try {
        H5XX_CALL(ret = H5Sclose(id));
        BOOST_FAIL("no h5xx::error exception thrown");
    }
    catch (h5xx::error const& e) {
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
        BOOST_REQUIRE(0 == strcmp(e.what(), "H5Sclose: not a data space"));
#else
        BOOST_REQUIRE(0 == strcmp(e.what(), "H5Sclose: not a dataspace"));
#endif
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
