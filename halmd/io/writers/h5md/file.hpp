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

#ifndef HALMD_IO_WRITERS_H5MD_FILE_HPP
#define HALMD_IO_WRITERS_H5MD_FILE_HPP

#include <boost/array.hpp>
#include <h5xx/h5xx.hpp>
#include <lua.hpp>
#include <string>

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

/**
 * H5MD file writer
 *
 * This class provides a common base for all H5MD file writers.
 * It creates the H5MD file and writes the H5MD metadata.
 */
class file
{
public:
    /** H5MD major and minor file version type */
    typedef boost::array<int, 2> version_type;

    /**
     * create H5MD file
     *
     * If author_name is an empty string it is retrieved from the password file
     * entry for the real user id of the calling process. If author_email is
     * empty output of this optional field is skipped.
     */
    file(std::string const& path, std::string const& author_name = "", std::string const& author_email = "");

    /** flush file to disk */
    void flush();
    /** explicitly close file */
    void close();
    /** get HDF5 root group */
    H5::Group root() const;
    /** get file pathname */
    std::string path() const;

    /** get H5MD file version */
    static version_type version();

    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    /** H5MD file */
    H5::H5File file_;
};

} // namespace h5md
} // namespace writers
} // namespace io
} // namespace halmd

#endif /* ! HALMD_IO_WRITERS_H5MD_FILE_HPP */
