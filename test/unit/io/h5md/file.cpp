/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2011 Peter Colberg
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

#define BOOST_TEST_MODULE file
#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>
#include <ctime>
#include <memory>

#include <halmd/io/writers/h5md/file.hpp>
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace halmd;
using namespace halmd::io; // avoid ambiguity of io:: between halmd::io and boost::io
using namespace std;

struct create_file
{
    typedef writers::h5md::file file_type;
    typedef file_type::version_type version_type;
    /** create H5MD file */
    create_file() {
        file = std::make_shared<file_type>("h5md.h5", "", "info@halmd.org");
    }
    /** close and unlink H5MD file */
    ~create_file() {
        file.reset();
        filesystem::remove("h5md.h5");
    }
    std::shared_ptr<file_type> file;
};

BOOST_FIXTURE_TEST_CASE( check_version, create_file )
{
    H5::Group group = file->root().openGroup("h5md");
    version_type version = h5xx::read_attribute<version_type>(group, "version");
    BOOST_CHECK_EQUAL( version[0], file->version()[0] );
    BOOST_CHECK_EQUAL( version[1], file->version()[1] );
}

BOOST_FIXTURE_TEST_CASE( read_attributes, create_file )
{
    H5::Group group = file->root().openGroup("h5md");
    version_type version = h5xx::read_attribute<version_type>(group, "version");
    BOOST_TEST_MESSAGE( "H5MD major version:\t" << version[0] );
    BOOST_TEST_MESSAGE( "H5MD minor version:\t" << version[1] );
    H5::Group creator = group.openGroup("creator");
    BOOST_TEST_MESSAGE( "H5MD creator name:\t\t" << h5xx::read_attribute<string>(creator, "name") );
    BOOST_TEST_MESSAGE( "H5MD creator version:\t" << h5xx::read_attribute<string>(creator, "version") );
    H5::Group author = group.openGroup("author");
    BOOST_TEST_MESSAGE( "H5MD author name:\t\t" << h5xx::read_attribute<string>(author, "name") );
    if (h5xx::exists_attribute(author, "email")) {
        BOOST_TEST_MESSAGE( "H5MD author email:\t\t" << h5xx::read_attribute<string>(author, "email") );
    }
}

BOOST_FIXTURE_TEST_CASE( flush_file, create_file )
{
    BOOST_CHECK_NO_THROW( file->flush() );
}

BOOST_FIXTURE_TEST_CASE( close_file, create_file )
{
    BOOST_CHECK_NO_THROW( file->close() );
}
