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
#include <ctime>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/writers/h5md/file.hpp>
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/realname.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

/**
 * This overwrites any existing file at the given path.
 */
file::file(string const& path)
{
    file_ = H5::H5File(path, H5F_ACC_TRUNC);

    H5::Group attr = file_.createGroup("h5md");
    h5xx::write_attribute(attr, "creation_time", time(NULL));
    h5xx::write_attribute(attr, "creator", PROGRAM_NAME);
    h5xx::write_attribute(attr, "creator_version", PROGRAM_VERSION);
    h5xx::write_attribute(attr, "version", file::version());
    h5xx::write_attribute(attr, "author", file::author());

    LOG("write to H5MD file: " << absolute_path(file_.getFileName()));
}

void file::flush()
{
    LOG("flush H5MD file: " << absolute_path(file_.getFileName()));
    file_.flush(H5F_SCOPE_GLOBAL);
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

file::version_type file::version()
{
    version_type version = {{ 0, 0 }};
    return version;
}

std::string file::author()
{
    return realname();
}

static signal<void ()>::slot_function_type
wrap_flush(shared_ptr<file> instance)
{
    return bind(&file::flush, instance);
}

void file::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("writers")
            [
                namespace_("h5md")
                [
                    class_<file, shared_ptr<file> >("file")
                        .def(constructor<string const&>())
                        // wrap as slot to support periodic flushing of file
                        .property("flush", &wrap_flush)
                        .def("close", &file::close)
                        .property("root", &file::root)
                        .property("path", &file::path)
                        .scope
                        [
                            def("version", &file::version)
                        ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_writers_h5md_file(lua_State* L)
{
    file::luaopen(L);
    return 0;
}

} // namespace h5md
} // namespace writers
} // namespace io
} // namespace halmd
