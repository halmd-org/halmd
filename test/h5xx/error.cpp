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

// #define H5XX_DEBUG
#include <h5xx/error.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(test_success)
{
    hid_t ret, id;
    H5XX_CHECK(ret = id = H5Screate(H5S_SCALAR));
    BOOST_REQUIRE(ret >= 0);
    H5XX_CHECK(ret = H5Sclose(id));
    BOOST_REQUIRE(ret >= 0);
}

BOOST_AUTO_TEST_CASE(test_throw_exception)
{
    hid_t ret, id;
    H5XX_CHECK(ret = id = H5Screate(H5S_SCALAR));
    H5XX_CHECK(ret = H5Sclose(id));
    try {
        H5XX_CHECK(ret = H5Sclose(id));
        BOOST_FAIL("no h5xx::error exception thrown");
    }
    catch (h5xx::error const&) {
        BOOST_REQUIRE(ret < 0);
    }
}

BOOST_AUTO_TEST_CASE(test_no_print)
{
    H5E_auto_t efunc = reinterpret_cast<H5E_auto_t>(H5Eprint);
    void* edata = stderr;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(H5Eget_auto(&efunc, &edata));
#else
    H5XX_CHECK(H5Eget_auto(H5E_DEFAULT, &efunc, &edata));
#endif
    BOOST_REQUIRE(NULL == efunc);
    BOOST_REQUIRE(NULL == edata);
}

BOOST_AUTO_TEST_CASE(test_reset_print)
{
    H5E_auto_t efunc = NULL;
#ifdef H5XX_USE_16_API
    void* edata = NULL;
    H5Eget_auto(&efunc, &edata);
#else
    void* edata = stderr;
    H5Eget_auto(H5E_DEFAULT, &efunc, &edata);
#endif
#ifdef H5_USE_16_API_DEFAULT
    BOOST_REQUIRE(reinterpret_cast<H5E_auto_t>(H5Eprint1) == efunc);
#else
    BOOST_REQUIRE(reinterpret_cast<H5E_auto_t>(H5Eprint) == efunc);
#endif
#ifdef H5XX_USE_16_API
    BOOST_REQUIRE(stderr == edata);
#else
    BOOST_REQUIRE(NULL == edata);
#endif
}

BOOST_AUTO_TEST_CASE(test_error_type)
{
    hid_t ret, id;
    H5XX_CHECK(ret = id = H5Screate(H5S_SCALAR));
    H5XX_CHECK(ret = H5Sclose(id));
    try {
        H5XX_CHECK(ret = H5Sclose(id));
        BOOST_FAIL("no h5xx::error exception thrown");
    }
    catch (h5xx::error const& e) {
#ifdef H5XX_USE_16_API
        BOOST_REQUIRE(e.stack.size() == 2);
        BOOST_CHECK(e.count(make_pair(H5E_ATOM, H5E_BADATOM)) == 1);
        BOOST_CHECK(e.stack[1] == make_pair(H5E_ATOM, H5E_BADATOM));
#else
        BOOST_REQUIRE(e.stack.size() == 1);
#endif
        BOOST_CHECK(e.count(make_pair(H5E_ARGS, H5E_BADTYPE)) == 1);
        BOOST_CHECK(e.stack[0] == make_pair(H5E_ARGS, H5E_BADTYPE));
    }
}

BOOST_AUTO_TEST_CASE(test_error_desc)
{
    hid_t ret, id;
    H5XX_CHECK(ret = id = H5Screate(H5S_SCALAR));
    H5XX_CHECK(ret = H5Sclose(id));
    try {
        H5XX_CHECK(ret = H5Sclose(id));
        BOOST_FAIL("no h5xx::error exception thrown");
    }
    catch (h5xx::error const& e) {
#ifdef H5XX_USE_16_API
        BOOST_REQUIRE_EQUAL(e.what(), "H5Sclose(): not a data space");
#else
        BOOST_REQUIRE_EQUAL(e.what(), "H5Sclose(): not a dataspace");
#endif
    }
}
