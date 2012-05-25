/*
 * Copyright © 2011  Peter Colberg
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

#include <boost/shared_ptr.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/readers/h5md/file.hpp>
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

bool file::check(string const& path)
{
    if (!H5::H5File::isHdf5(path)) {
        return false;
    }
    H5::H5File file(path, H5F_ACC_RDONLY);
    if (!h5xx::exists_group(file, "h5md")) {
        return false;
    }
    return true;
}

file::file(string const& path)
{
    file_ = H5::H5File(path, H5F_ACC_RDONLY);

    H5::Group attr = file_.openGroup("h5md");
    creation_time_ = h5xx::read_attribute<time_t>(attr, "creation_time");
    creator_ = h5xx::read_attribute<string>(attr, "creator");
    creator_version_ = h5xx::read_attribute<string>(attr, "creator_version");
    version_ = h5xx::read_attribute<version_type>(attr, "version");
    author_ = h5xx::read_attribute<string>(attr, "author");

    LOG("read from H5MD file: " << absolute_path(file_.getFileName()));
}

void file::close()
{
    file_.close();
}

H5::Group file::root() const
{
    return file_.openGroup("/");
}

string file::path() const
{
    return file_.getFileName();
}

void file::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("readers")
            [
                namespace_("h5md")
                [
                    class_<file, boost::shared_ptr<file> >("file")
                        .def(constructor<string const&>())
                        .def("close", &file::close)
                        .property("root", &file::root)
                        .property("path", &file::path)
                        .property("version", &file::version)
                        .property("creator", &file::creator)
                        .property("creator_version", &file::creator_version)
                        .property("creation_time", &file::creation_time)
                        .property("author", &file::author)
                        .scope
                        [
                            def("check", &file::check)
                        ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_readers_h5md_file(lua_State* L)
{
    file::luaopen(L);
    return 0;
}

} // namespace h5md
} // namespace readers
} // namespace io
} // namespace halmd
