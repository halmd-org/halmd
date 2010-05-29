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
#define BOOST_TEST_MODULE test_h5xx_file
#include <boost/test/unit_test.hpp>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <fstream>

// #define H5XX_DEBUG
#include <h5xx/file.hpp>

#define H5XX_TEST_FILENAME "test_h5xx_file.h5"

using namespace boost;
using namespace boost::filesystem;
using namespace std;

#define BOOST_CHECK_H5XX_EXCEPTION(expression, major, minor) \
    BOOST_CHECK_EXCEPTION( \
        expression \
      , h5xx::error \
      , bind(&h5xx::error::count, _1, make_pair(major, minor)) \
    )

struct remove_file
{
    // remove test file before and after test case
    remove_file()  { remove(H5XX_TEST_FILENAME); }
    ~remove_file() { remove(H5XX_TEST_FILENAME); }
};

BOOST_FIXTURE_TEST_CASE(file_format_is_hdf5, remove_file)
{
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file::format(H5XX_TEST_FILENAME), H5E_FILE, H5E_CANTOPENFILE);
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file(H5XX_TEST_FILENAME), H5E_FILE, H5E_CANTOPENFILE);
    (void) h5xx::file(H5XX_TEST_FILENAME, h5xx::file::excl);
    BOOST_CHECK(h5xx::file::format(H5XX_TEST_FILENAME));
}

BOOST_FIXTURE_TEST_CASE(file_format_is_not_hdf5, remove_file)
{
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file::format(H5XX_TEST_FILENAME), H5E_FILE, H5E_CANTOPENFILE);
    ofstream file;
    file.exceptions(ofstream::eofbit | ofstream::failbit | ofstream::badbit);
    file.open(H5XX_TEST_FILENAME);
    file.close();
    BOOST_CHECK(!h5xx::file::format(H5XX_TEST_FILENAME));
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file(H5XX_TEST_FILENAME), H5E_FILE, H5E_NOTHDF5);
    file.open(H5XX_TEST_FILENAME);
    file << "Hello, World!" << std::endl;
    file.close();
    BOOST_CHECK(!h5xx::file::format(H5XX_TEST_FILENAME));
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file(H5XX_TEST_FILENAME), H5E_FILE, H5E_NOTHDF5);
}

BOOST_FIXTURE_TEST_CASE(file_open_modes, remove_file)
{
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::ro), H5E_FILE, H5E_CANTOPENFILE);
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::rw), H5E_FILE, H5E_CANTOPENFILE);
    BOOST_CHECK_NO_THROW(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::excl));
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::excl), H5E_FILE, H5E_CANTOPENFILE);
    BOOST_CHECK_NO_THROW(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::trunc));
    BOOST_CHECK_NO_THROW(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::ro));
    BOOST_CHECK_NO_THROW(h5xx::file(H5XX_TEST_FILENAME, h5xx::file::rw));
}

BOOST_FIXTURE_TEST_CASE(file_name, remove_file)
{
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    BOOST_CHECK_EQUAL(file.name(), H5XX_TEST_FILENAME);
}

BOOST_FIXTURE_TEST_CASE(file_flush, remove_file)
{
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    BOOST_CHECK_NO_THROW(file.flush());
}

BOOST_FIXTURE_TEST_CASE(file_open_read_only, remove_file)
{
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    hid_t hid1;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(hid1 = H5Gcreate(file.hid(), "hello", 0));
#else
    H5XX_CHECK(hid1 = H5Gcreate(file.hid(), "hello", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
#endif
    H5XX_CHECK(H5Gclose(hid1));

    file = h5xx::file(H5XX_TEST_FILENAME);
    hid_t hid2;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(hid2 = H5Gopen(file.hid(), "hello"));
#else
    H5XX_CHECK(hid2 = H5Gopen(file.hid(), "hello", H5P_DEFAULT));
#endif
    H5XX_CHECK(H5Gclose(hid2));
    hid_t hid3;
#ifdef H5XX_USE_16_API
    BOOST_CHECK_H5XX_EXCEPTION(H5XX_CHECK(hid3 = H5Gcreate(file.hid(), "world", 0)), H5E_SYM, H5E_CANTINIT);
#else
    BOOST_CHECK_H5XX_EXCEPTION(H5XX_CHECK(hid3 = H5Gcreate(file.hid(), "world", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)), H5E_SYM, H5E_CANTINIT);
#endif
}

BOOST_FIXTURE_TEST_CASE(file_open_read_write, remove_file)
{
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    hid_t hid1;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(hid1 = H5Gcreate(file.hid(), "hello", 0));
#else
    H5XX_CHECK(hid1 = H5Gcreate(file.hid(), "hello", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
#endif
    H5XX_CHECK(H5Gclose(hid1));

    file = h5xx::file(H5XX_TEST_FILENAME, h5xx::file::rw);
    hid_t hid2;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(hid2 = H5Gopen(file.hid(), "hello"));
#else
    H5XX_CHECK(hid2 = H5Gopen(file.hid(), "hello", H5P_DEFAULT));
#endif
    H5XX_CHECK(H5Gclose(hid2));
    hid_t hid3;
#ifdef H5XX_USE_16_API
    BOOST_CHECK_NO_THROW(H5XX_CHECK(hid3 = H5Gcreate(file.hid(), "world", 0)));
#else
    BOOST_CHECK_NO_THROW(H5XX_CHECK(hid3 = H5Gcreate(file.hid(), "world", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)));
#endif
}
