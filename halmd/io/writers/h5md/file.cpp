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

#include <ctime>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/writers/h5md/file.hpp>
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/realname.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/version.h>

using namespace std;

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

file::file(string const& path, string const& author_name, string const& author_email, bool overwrite)
{
    if (boost::filesystem::exists(path)) {
        if (overwrite) {
            LOG_WARNING("overwrite existing file \"" << path << "\"");
        } else {
            throw std::runtime_error("Output file \"" + path + "\" exists and won't be overwritten");
        }
    }

    // open file with write access, truncate file if it exists
    file_ = H5::H5File(path, H5F_ACC_TRUNC);

    H5::Group h5md = file_.createGroup("h5md");
    h5xx::write_attribute(h5md, "version", file::version());

    H5::Group creator = h5md.createGroup("creator");
    h5xx::write_attribute(creator, "name", PROGRAM_DESC_ASCII);
    h5xx::write_attribute(creator, "version", PROGRAM_VERSION PROGRAM_VARIANT);

    H5::Group author = h5md.createGroup("author");
    h5xx::write_attribute(author, "name", author_name.empty() ? realname() : author_name);
    if (!author_email.empty()) {
        h5xx::write_attribute(author, "email", author_email);
    }

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
    version_type version = {{ 1, 0 }};
    return version;
}

void file::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("writers")
            [
                namespace_("h5md")
                [
                    class_<file, std::shared_ptr<file> >("file")
                        .def(constructor<string const&, string const&, string const&, bool>())
                        .def("flush", &file::flush)
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
