/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2011 Peter Colberg
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

#ifndef HALMD_IO_READERS_H5MD_FILE_HPP
#define HALMD_IO_READERS_H5MD_FILE_HPP

#include <boost/array.hpp>
#include <h5xx/h5xx.hpp>
#include <lua.hpp>
#include <string>

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

/**
 * H5MD file reader
 *
 * This class opens the H5MD file and reads the H5MD metadata.
 */
class file
{
public:
    /** H5MD major and minor file version type */
    typedef boost::array<int, 2> version_type;

    /** check that file is H5MD file */
    static bool check(std::string const& path);
    /** open H5MD file in read-only mode */
    file(std::string const& path);
    /** explicitly close file */
    void close();
    /** get HDF5 root group */
    H5::Group root() const;
    /** get file pathname */
    std::string path() const;
    /** Lua bindings */
    static void luaopen(lua_State* L);

    /** get H5MD file version */
    version_type const& version() const
    {
        return version_;
    }

    /** get name of H5MD file creator */
    std::string const& creator_name() const
    {
        return creator_name_;
    }

    /** get version of H5MD file creator */
    std::string const& creator_version() const
    {
        return creator_version_;
    }

    /** get name of H5MD file author */
    std::string const& author_name() const
    {
        return author_name_;
    }

    /** get name of H5MD file author */
    std::string const& author_email() const
    {
        return author_email_;
    }

private:
    /** H5MD file */
    H5::H5File file_;
    /** H5MD file version */
    version_type version_;
    /** H5MD meta data */
    std::string creator_name_;
    std::string creator_version_;
    std::string author_name_;
    std::string author_email_;
};

} // namespace h5md
} // namespace readers
} // namespace io
} // namespace halmd

#endif /* ! HALMD_IO_READERS_H5MD_FILE_HPP */
