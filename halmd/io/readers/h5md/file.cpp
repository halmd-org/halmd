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

#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/readers/h5md/file.hpp>
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/version.h>

using namespace std;

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

template <typename T>
bool exists_attribute_of_type(H5::Group const& group, string name, bool mandatory = true)
{
    if (h5xx::exists_attribute(group, name)) {
        return h5xx::has_type<T>(group.openAttribute(name));
    }
    else {
        return !mandatory;
    }
}

bool file::check(string const& path)
{
    if (!H5::H5File::isHdf5(path)) {
        return false;
    }
    H5::H5File file(path, H5F_ACC_RDONLY);

    if (!h5xx::exists_group(file, "h5md")) {
        return false;
    }
    H5::Group h5md = file.openGroup("h5md");

    if(!exists_attribute_of_type<version_type>(h5md, "version")) {
        return false;
    }

    if (h5xx::exists_group(h5md, "creator")) {
        H5::Group group = h5md.openGroup("creator");
        if(!exists_attribute_of_type<string>(group, "name") ||
           !exists_attribute_of_type<string>(group, "version")) {
            return false;
        }
    }
    else {
        return false;
    }

    if (h5xx::exists_group(h5md, "author")) {
        H5::Group group = h5md.openGroup("author");
        if(!exists_attribute_of_type<string>(group, "name") ||
           !exists_attribute_of_type<string>(group, "email", false)) {
            return false;
        }
    }
    else {
        return false;
    }

    return true;
}

file::file(string const& path)
{
    file_ = H5::H5File(path, H5F_ACC_RDONLY);

    H5::Group h5md = file_.openGroup("h5md");
    version_ = h5xx::read_attribute<version_type>(h5md, "version");

    H5::Group creator = h5md.openGroup("creator");
    creator_name_ = h5xx::read_attribute<string>(creator, "name");
    creator_version_ = h5xx::read_attribute<string>(creator, "version");

    H5::Group author = h5md.openGroup("author");
    author_name_ = h5xx::read_attribute<string>(author, "name");
    if (h5xx::exists_attribute(author, "email")) {
        author_email_ = h5xx::read_attribute<string>(author, "email");
    }

    LOG("read from H5MD file: " << absolute_path(file_.getFileName()));

    if (version_ != version_type({{ 1, 0 }})) {
        LOG_ERROR("version of H5MD format: " << version_[0] << "." << version_[1]);
        throw std::logic_error("unsupported H5MD file");
    }
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
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("readers")
            [
                namespace_("h5md")
                [
                    class_<file, std::shared_ptr<file> >("file")
                        .def(constructor<string const&>())
                        .def("close", &file::close)
                        .property("root", &file::root)
                        .property("path", &file::path)
                        .property("version", &file::version)
                        .property("creator_name", &file::creator_name)
                        .property("creator_version", &file::creator_version)
                        .property("author_name", &file::author_name)
                        .property("author_email", &file::author_email)
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
