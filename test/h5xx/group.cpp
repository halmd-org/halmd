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

#include <boost/foreach.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_h5xx_group
#include <boost/test/unit_test.hpp>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/test/output_test_stream.hpp>

// #define H5XX_DEBUG
#include <h5xx/file.hpp>
#include <h5xx/group.hpp>

#define H5XX_TEST_FILENAME "test_h5xx_group.h5"

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

BOOST_FIXTURE_TEST_CASE(test_group_operator, remove_file)
{
    boost::test_tools::output_test_stream os;
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    os << file;
    BOOST_CHECK(os.is_equal("/"));
    os << file / "hello";
    BOOST_CHECK(os.is_equal("/hello"));
    os << file / "hello" / "," / "world" / "!";
    BOOST_CHECK(os.is_equal("/hello/,/world/!"));
    os << file / "another" / "hello" / "," / "world" / "!";
    BOOST_CHECK(os.is_equal("/another/hello/,/world/!"));
    os << file / "another" / "hello" / "," / "world" / "!" / "././././:" / ";";
    BOOST_CHECK(os.is_equal("/another/hello/,/world/!/:/;"));
    os << file / "a n o t h e r" / "h e l l o" / "," / " W O R L D " / "!" / "/:" / ";";
    BOOST_CHECK(os.is_equal("/:/;"));
    BOOST_CHECK(h5xx::group::exists(file / "a n o t h e r" / "h e l l o" / "," / " W O R L D ", "!"));
    BOOST_CHECK(h5xx::group::exists(file / "a n o t h e r/h e l l o/,/ W O R L D ", "!"));
    BOOST_CHECK(h5xx::group::exists(file / "a n o t h e r/h e l l o/,/ W O R L D ", "."));
    BOOST_CHECK(!h5xx::group::exists(file / "a n o t h e r/h e l l o/,/ W O R L D ", "-"));
    BOOST_CHECK(h5xx::group::exists(file / ":", ";"));
    BOOST_CHECK(!h5xx::group::exists(file / ":", "-"));

    file = h5xx::file(H5XX_TEST_FILENAME);
    os << file;
    BOOST_CHECK(os.is_equal("/"));
    os << file / "hello";
    BOOST_CHECK(os.is_equal("/hello"));
    os << file / "hello" / ",";
    BOOST_CHECK(os.is_equal("/hello/,"));
    os << file / "hello" / "," / "world";
    BOOST_CHECK(os.is_equal("/hello/,/world"));
    os << file / "hello" / "," / "world" / "!";
    BOOST_CHECK(os.is_equal("/hello/,/world/!"));
    os << file / "another" / "hello" / "," / "world" / "!";
    BOOST_CHECK(os.is_equal("/another/hello/,/world/!"));
    os << file / "another" / "hello" / "," / "world" / "!" / "././././:" / ";";
    BOOST_CHECK(os.is_equal("/another/hello/,/world/!/:/;"));
    os << file / "a n o t h e r" / "h e l l o" / "," / " W O R L D " / "!";
    BOOST_CHECK(os.is_equal("/a n o t h e r/h e l l o/,/ W O R L D /!"));
    os << file / "a n o t h e r" / "h e l l o" / "," / " W O R L D " / "!" / "/:" / ";";
    BOOST_CHECK(os.is_equal("/:/;"));
    BOOST_CHECK_H5XX_EXCEPTION(file / "does", H5E_SYM, H5E_CANTINIT);
    BOOST_CHECK_H5XX_EXCEPTION(file / "does", H5E_SYM, H5E_CANTINIT);
    BOOST_CHECK_H5XX_EXCEPTION(file / "does" / "not" / "exist", H5E_SYM, H5E_CANTINIT);
    BOOST_CHECK_H5XX_EXCEPTION(file / "does" / "not" / "exist", H5E_SYM, H5E_CANTINIT);
}

BOOST_FIXTURE_TEST_CASE(test_group_assignment, remove_file)
{
    boost::test_tools::output_test_stream os;
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    BOOST_CHECK_EQUAL(file.path(), "/");
    h5xx::group group = file / "/";
    BOOST_CHECK_EQUAL(group.path(), "/");
    string path = "";
    for (char i = 'a'; i <= 'z'; ++i) {
        group /= string(1, i);
        BOOST_CHECK_EQUAL(group.path(), path += "/" + string(1, i));
    }
    BOOST_CHECK_EQUAL(group.path(), "/a/b/c/d/e/f/g/h/i/j/k/l/m/n/o/p/q/r/s/t/u/v/w/x/y/z");
    group = group / "/" / "_";
    path = "/_";
    for (char i = 'A'; i <= 'Z'; ++i) {
        group /= string(1, i);
        BOOST_CHECK_EQUAL(group.path(), path += "/" + string(1, i));
    }
    BOOST_CHECK_EQUAL(group.path(), "/_/A/B/C/D/E/F/G/H/I/J/K/L/M/N/O/P/Q/R/S/T/U/V/W/X/Y/Z");
}

BOOST_FIXTURE_TEST_CASE(test_group_iterator, remove_file)
{
    boost::test_tools::output_test_stream os;
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    for (char i = 'A'; i <= 'Z'; ++i) {
        file / string(1, i);
    }
    h5xx::group::iterator const end;
    size_t children = 0;
    for (h5xx::group::iterator it = file; it != end; ++children) {
        os << (*it / "visited");
        BOOST_CHECK(os.is_equal((*it).path() + "/visited"));
        BOOST_CHECK(it != it++);
    }
    BOOST_CHECK_EQUAL(children, 26u);
    for (char i = 'A'; i <= 'Z'; ++i) {
        file / string(1, i);
    }

    file = h5xx::file(H5XX_TEST_FILENAME);
    for (char i = 'A'; i <= 'Z'; ++i) {
        BOOST_CHECK(h5xx::group::exists(file / string(1, i), "visited"));
    }
}

BOOST_FIXTURE_TEST_CASE(test_group_iterator_special_cases, remove_file)
{
    boost::test_tools::output_test_stream os;
    h5xx::file file(H5XX_TEST_FILENAME, h5xx::file::excl);
    size_t children = 0;
    h5xx::group::iterator /* omit const */ end;
    for (h5xx::group::iterator it = file; it != end; ++it, ++children) {
        os << (*it / "visited");
        BOOST_CHECK(os.is_equal((*it).path() + "/visited"));
    }
    BOOST_CHECK_EQUAL(children, 0u);
    hid_t hid1;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(hid1 = H5Dcreate(file.hid(), "c", H5T_NATIVE_FLOAT, H5Screate(H5S_SCALAR), H5P_DEFAULT));
#else
    H5XX_CHECK(hid1 = H5Dcreate(file.hid(), "c", H5T_NATIVE_FLOAT, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
#endif
    H5XX_CHECK(H5Dclose(hid1));
#ifdef H5XX_USE_16_API
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::group::exists(file, "c"), H5E_ARGS, H5E_BADTYPE);
#else
    BOOST_CHECK_H5XX_EXCEPTION(h5xx::group::exists(file, "c"), H5E_SYM, H5E_BADTYPE);
#endif
    children = 0;
    for (h5xx::group::iterator it = file; it != end; ++it, ++children) {
        os << (*it / "visited");
        BOOST_CHECK(os.is_equal((*it).path() + "/visited"));
    }
    BOOST_CHECK_EQUAL(children, 0u);
    (void) (file / "b");
    children = 0;
    for (h5xx::group::iterator it = file; it != end; ++it, ++children) {
        os << (*it / "visited");
        BOOST_CHECK(os.is_equal("/b/visited"));
    }
    BOOST_CHECK_EQUAL(children, 1u);
    hid_t hid2;
#ifdef H5XX_USE_16_API
    H5XX_CHECK(hid2 = H5Dcreate(file.hid(), "a", H5T_NATIVE_FLOAT, H5Screate(H5S_SCALAR), H5P_DEFAULT));
#else
    H5XX_CHECK(hid2 = H5Dcreate(file.hid(), "a", H5T_NATIVE_FLOAT, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
#endif
    H5XX_CHECK(H5Dclose(hid2));
    children = 0;
    for (h5xx::group::iterator it = file; it != end; ++it, ++children) {
        BOOST_CHECK_EQUAL(it->path(), "/b");
    }
    BOOST_CHECK_EQUAL(children, 1u);
}
